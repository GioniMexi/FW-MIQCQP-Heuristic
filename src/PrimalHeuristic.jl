module PrimalHeuristic

using ArgParse
using Boscia
using Dates
using FrankWolfe
using Graphs
using Gurobi
using JSON
using LinearAlgebra
import MathOptInterface as MOI
using SparseArrays
using ProgressMeter
using TimerOutputs

include("presolvers.jl")
include("simple_stats.jl")
include("simple_parser.jl")
include("best_known_sol.jl")
include("propagation.jl")
include("fix_and_prop.jl")
include("lns_heuristics.jl")
include("variable_bounds.jl")
include("utils.jl")
include("main.jl")
include("HeuristicsStructs.jl")
include("post_heuristics_callback.jl")
include("penalty_functions.jl")
include("parallelisation_strategy.jl")

"""
Construct Boscia arguments: modified objective function, gradient, and heuristics.
"""
function build_objective_and_heuristics(problem_data::ProblemData, heur_params::HeurParams, lmo)
    Q, b, c0 = problem_data.quad_objective_data

    quad_obj = !isempty(Q)
    quad_cons = length(problem_data.quadratic_lowerthan) + length(problem_data.quadratic_equalto) > 0

    convexify = false
    shift = 0.0
    p = is_type_qbo(problem_data.type) ? heur_params.convexify_objective : 0.0
    if p > 0
        # Convexifying with parameter p
        relevant_index = min(problem_data.heur_stats.num_vars, Int(floor(problem_data.heur_stats.num_vars * (1 - p))) + 1)
        if problem_data.spectrum[relevant_index] < 0
            convexify = true
            shift = -problem_data.spectrum[relevant_index]
        end
    end

    # Define the objective function
    function f(x)
        obj = 0.0
        # quadratic objective
        if quad_obj
            obj += 0.5 * dot(x, Q, x) + dot(b, x) + c0
            # linear objective
        else
            obj += dot(b, x) + c0
        end
        if convexify
            obj += 0.5 * shift * (dot(x, x) - sum(x))
        end
        # quadratic constraints
        if quad_cons
            if heur_params.power_penalty > 0
                obj += objective_power_penalty(x, problem_data, heur_params.mu_barrier, heur_params.power_penalty)
            else
                obj += objective_log_barrier_penalty(x, problem_data, heur_params.mu_barrier, heur_params.shift_barrier)
            end
        end
        return obj
    end

    # Define the gradient of the objective function
    function grad!(storage, x)
        if quad_obj
            mul!(storage, Q, x)  # add Qx to the `storage`
            storage .+= b        # add the linear coefficients vector b to the `storage`
        else
            copyto!(storage, b) # copy the coefficients vector b to the `storage`
        end
        if convexify
            @. storage += shift * (x - 0.5)
        end
        if quad_cons
            if heur_params.power_penalty > 0
                gradient_power_penalty!(storage, x, problem_data, heur_params.mu_barrier, heur_params.power_penalty)
            else
                gradient_log_barrier_penalty!(storage, x, problem_data, heur_params.mu_barrier, heur_params.shift_barrier)
            end
        end
        return nothing
    end

    # Create the post heuristic callback
    post_heuristics_callback = build_post_heuristic_callback(problem_data, heur_params)

    line_search = heur_params.use_secant ? FrankWolfe.Secant() : FrankWolfe.Adaptive()

    # Create the RENS heuristic
    rens_heu = Boscia.Heuristic((tree, blmo, x) -> rens_heuristic(tree, blmo, x, problem_data, heur_params), heur_params.prob_rens, :rens)
    # Create the Active Set RENS heuristic
    active_set_rens_heu = Boscia.Heuristic((tree, blmo, x) -> active_set_rens_heuristic(tree, blmo, x, problem_data, heur_params), heur_params.prob_arens, :asrens)
    # Create the RINS heuristic
    rins_heu = Boscia.Heuristic((tree, blmo, x) -> rins_heuristic(tree, blmo, x, problem_data, heur_params), heur_params.prob_rins, :rins)
    # Create Undercover Heuristic
    undercover_heu = Boscia.Heuristic((tree, blmo, x) -> undercover_heuristic(tree, blmo, x, problem_data, heur_params), heur_params.prob_undercover, :undercover)
    # Create Alternating Heuristic
    alternating_heu = Boscia.Heuristic((tree, blmo, x) -> alternating_heuristic(tree, blmo, x, problem_data, heur_params), heur_params.prob_alternating, :alternating)

    hyperplane_aware_rounding = typeof(lmo) in [Boscia.UnitSimplexSimpleBLMO, Boscia.ProbabilitySimplexSimpleBLMO, Boscia.ReverseKnapsackBLMO] ? heur_params.prob_custom_heu : 0.0
    hyperplane_aware_heu = Boscia.Heuristic(Boscia.rounding_hyperplane_heuristic, hyperplane_aware_rounding, :hyperplane_aware_rounding)

    follow_grad_heu = Boscia.Heuristic((tree, blmo, x) -> Boscia.follow_gradient_heuristic(tree, blmo, x, heur_params.max_iter_follow_gradient_heu), heur_params.prob_follow_gradient_heu, :follow_grad)

    rounding_01_prob = problem_data.heur_stats.num_int_vars == 0 ? heur_params.prob_rounding_01_heu : 0.0
    rounding_01_heu = Boscia.Heuristic(Boscia.rounding_lmo_01_heuristic, rounding_01_prob, :rounding01)

    probability_rounding_prob = problem_data.heur_stats.num_int_vars == 0 ? heur_params.prob_probability_rounding_heu : 0.0
    probability_rounding_heu = Boscia.Heuristic(Boscia.probability_rounding, probability_rounding_prob, :prob_rounding)

    heuristics_for_Boscia = [rens_heu, active_set_rens_heu, rins_heu, undercover_heu, alternating_heu, hyperplane_aware_heu, follow_grad_heu, rounding_01_heu, probability_rounding_heu]

    return f, grad!, post_heuristics_callback, heuristics_for_Boscia, line_search
end

"""
Run the heuristic.
"""
function run_heuristic(problem_data::ProblemData, heur_params::HeurParams, optimizer_linear, i, nthreads)
    num_vars = problem_data.heur_stats.num_vars
    N = problem_data.N
    # define lmo
    if problem_data.type in [QUBO, QUBO_BIPARTITE]
        heur_params.prob_custom_heu = heur_params.prob_custom_heu == 0.0 ? 1.0 : heur_params.prob_custom_heu
        lmo = Boscia.ManagedBoundedLMO(Boscia.CubeSimpleBLMO(fill(0.0, num_vars), fill(1.0, num_vars), collect(1:num_vars)),
            fill(0.0, num_vars),
            fill(1.0, num_vars),
            collect(1:num_vars),
            num_vars
        )
    elseif problem_data.type == QBO_KNAPSACK_UNIT
        heur_params.prob_custom_heu = heur_params.prob_custom_heu == 0.0 ? 1.0 : heur_params.prob_custom_heu
        lmo = Boscia.ManagedBoundedLMO(Boscia.UnitSimplexSimpleBLMO(N),
            fill(0.0, num_vars),
            fill(1.0, num_vars),
            collect(1:num_vars),
            num_vars
        )
    elseif problem_data.type == QBO_KNAPSACK_PROB
        heur_params.prob_custom_heu = heur_params.prob_custom_heu == 0.0 ? 1.0 : heur_params.prob_custom_heu
        lmo = Boscia.ManagedBoundedLMO(Boscia.ProbabilitySimplexSimpleBLMO(N),
            fill(0.0, num_vars),
            fill(1.0, num_vars),
            collect(1:num_vars),
            num_vars
        )
    elseif problem_data.type == QBO_KNAPSACK_UNIT_REV
        heur_params.prob_custom_heu = heur_params.prob_custom_heu == 0.0 ? 1.0 : heur_params.prob_custom_heu
        lmo = Boscia.ManagedBoundedLMO(Boscia.ReverseKnapsackBLMO(num_vars; N = N),
            fill(0.0, num_vars),
            fill(1.0, num_vars),
            collect(1:num_vars),
            num_vars
        )
    else
        lmo = FrankWolfe.MathOptLMO(optimizer_linear)
    end

    time_ref = problem_data.time_ref

    f, grad!, post_heuristics_callback, heuristics_for_Boscia, line_search = build_objective_and_heuristics(problem_data, heur_params, lmo)

    lmo_time_estimate = Inf
    if heur_params.perform_benchmark_oracle
        to = benchmark_oracles(f, grad!, () -> randn(problem_data.heur_stats.num_vars), lmo; k = 1, nocache = true)
        lmo_time_estimate = TimerOutputs.time(to["lmo"]) / 1e9
    end

    heur_params.fw_max_iter = isfinite(lmo_time_estimate) ? min(Int(floor(10 / (lmo_time_estimate))), 10000) : heur_params.fw_max_iter

    
    boscia_time_limit = round(Float64, remaining_time_limit(time_ref, heur_params.time_limit)) - 2.0



    # non-convex settings for Boscia
    no_pruning = true
    ignore_lower_bound = true
    add_all_solutions = true
    if heur_params.convexify_objective == 1.0
        no_pruning = false
        ignore_lower_bound = false
        add_all_solutions = false
    end

    if heur_params.fw_variant == 0
        fw_variant = Boscia.BPCG()
    elseif heur_params.fw_variant == 1
        fw_variant = Boscia.VanillaFrankWolfe()
    else
        @error "Unknown Frank-Wolfe variant"
    end

    x, tlmo, result = Boscia.solve(
        f, grad!, lmo; verbose = heur_params.verbose_boscia, fw_verbose = heur_params.verbose_fw, print_iter = 1,
        use_postsolve = false, no_pruning = no_pruning, ignore_lower_bound = ignore_lower_bound, add_all_solutions = add_all_solutions,
        post_heuristics_callback = post_heuristics_callback, custom_heuristics = heuristics_for_Boscia,
        time_limit = boscia_time_limit, node_limit = heur_params.node_limit, line_search = line_search, lazy = heur_params.fw_lazy,
        variant = fw_variant, max_fw_iter = heur_params.fw_max_iter, fw_timeout = heur_params.fw_timeout,
        max_time_lmo = heur_params.max_time_lmo
    )

    cnt = 0
    if heur_params.restart
        while true
            cnt += 1

            boscia_time = result[:total_time_in_sec]
            if remaining_time_limit(problem_data.time_ref, heur_params.time_limit) < min(heur_params.restart_min_time, boscia_time)
                break
            end
            prev_solutions = result[:tree_solutions]
            sort!(prev_solutions; by = s -> s.objective)

            previous_sol = !isempty(prev_solutions) ? prev_solutions[1].solution : nothing

            restart_strategy!(heur_params, problem_data, i, nthreads, cnt)

            f, grad!, post_heuristics_callback, heuristics_for_Boscia, line_search = build_objective_and_heuristics(problem_data, heur_params, lmo)

            gradient = similar(x)
            grad!(gradient, x)
            v = mod(cnt, 2) == 0 ? compute_extreme_point(tlmo.blmo, randn(length(x))) : compute_extreme_point(tlmo.blmo, gradient)
            active_set = FrankWolfe.ActiveSet([(1.0, v)])

            time_ref = problem_data.time_ref
            boscia_time_limit = round(Float64, remaining_time_limit(time_ref, heur_params.time_limit)) - 2.0
            
            x, tlmo, result = Boscia.solve(
                f, grad!, tlmo.blmo; verbose = heur_params.verbose_boscia, fw_verbose = heur_params.verbose_fw, print_iter = 1,
                use_postsolve = false, no_pruning = no_pruning, ignore_lower_bound = ignore_lower_bound, add_all_solutions = add_all_solutions,
                post_heuristics_callback = post_heuristics_callback, custom_heuristics = heuristics_for_Boscia,
                time_limit = boscia_time_limit, node_limit = heur_params.node_limit,
                line_search = line_search, lazy = heur_params.fw_lazy, active_set = active_set, start_solution = previous_sol,
                max_fw_iter = heur_params.fw_max_iter, fw_timeout = heur_params.fw_timeout, fw_variant = fw_variant,
                max_time_lmo = heur_params.max_time_lmo
            )
        end
    end
    @info "Primal heuristic finished on thread $i with $cnt restarts"
end

function main(; time_ref = Dates.now(), file_path, result_path = "", heur_params = HeurParams())
    @info "Starting primal heuristic"

    # Disable propagation if FixAndProp is not compiled
    if !FIXANDPROP_AVAILABLE
        if heur_params.use_propagation_global || heur_params.use_propagation_undercover
            @warn "FixAndProp not compiled. Disabling propagation (use_propagation_global=false, use_propagation_undercover=false)."
            heur_params.use_propagation_global = false
            heur_params.use_propagation_undercover = false
        end
    end

    # get file name without extension and path
    name = split(split(file_path, '/')[end], '.')[1]

    # Create the model
    original_model = create_moi_model(file_path)

    # Get the linear model and extract the quadratic constraints and the objective
    linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data = copy_mip_new_model(original_model, MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()))

    walltime_start = elapsed_time(time_ref)

    # Presolve the bilinear constraints (modifies first arguments in place)
    quad_objective_data, postsolve_bilinear = presolve_bilinear!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model)
    linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model, postsolve_deperspective =
        presolve_deperspective(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model; use_indicator = true)

    if heur_params.use_mccormick
        mccormick_original_model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
        mccormick_linear_model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
        MOI.copy_to(mccormick_original_model, original_model)
        MOI.copy_to(mccormick_linear_model, linear_model)

        num_added_mccormick_vars, num_added_mccormick_cons, quadratic_lowerthan_mccormick, quadratic_equalto_mccormick, quad_objective_data_mccormick, postsolve_mccormick = mccormick_relaxation!(mccormick_linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, mccormick_original_model, heur_params.bigM_McCormick)
        heur_params.use_mccormick = num_added_mccormick_cons > 0
    end

    postsolve_function = x -> postsolve_bilinear(postsolve_deperspective(x))

    sense = MOI.get(original_model, MOI.ObjectiveSense())
    initial_obj_value = sense == MOI.MIN_SENSE ? Inf : -Inf

    type = UNKNOWN

    num_vars = num_variables(linear_model)
    num_bin = length(MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}()))

    variable_lower_bounds = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{Float64}}())
    variable_upper_bounds = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{Float64}}())
    variable_equalities = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.EqualTo{Float64}}())
    bounds_present = !isempty(variable_lower_bounds) || !isempty(variable_upper_bounds) || !isempty(variable_equalities)

    less_than_constraints = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}())
    greater_than_constraints = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}}())
    equal_to_constraints = MOI.get(linear_model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}())
    constraints_present = !isempty(less_than_constraints) || !isempty(greater_than_constraints) || !isempty(equal_to_constraints)

    ind1 = Int64[]
    ind2 = Int64[]
    N = 0
    spectrum = Float64[]

    if num_vars == num_bin # binary problem
        if isempty(quad_objective_data[1]) # linear objective
            type = LBO

        else
            spectrum = eigvals(collect(quad_objective_data[1]))
            if !bounds_present && !constraints_present # unconstrained
                g = SimpleGraph(num_vars)
                for (i, j, v) in zip(findnz(quad_objective_data[1])...)
                    if i != j
                        add_edge!(g, i, j)
                    end
                end
                m = bipartite_map(g)
                if length(m) != num_vars
                    type = QUBO
                else
                    type = QUBO_BIPARTITE
                    heur_params.prob_alternating = heur_params.prob_alternating == 0.0 ? 1.0 : heur_params.prob_alternating
                    ind1 = findall(m .== 1)
                    ind2 = findall(m .== 2)
                end

            else
                type = QBO
                if isempty(quadratic_lowerthan) && isempty(quadratic_equalto) # no quadratic constraints
                    knapsack_constraint = false
                    s = :lessthan
                    if length(less_than_constraints) + length(greater_than_constraints) + length(equal_to_constraints) == 1
                        constraint, s = if !isempty(less_than_constraints)
                            (less_than_constraints[1], :lessthan)
                        elseif !isempty(greater_than_constraints)
                            (greater_than_constraints[1], :greaterthan)
                        else
                            (equal_to_constraints[1], :equalto)
                        end
                        func = MOI.get(linear_model, MOI.ConstraintFunction(), constraint)
                        set = MOI.get(linear_model, MOI.ConstraintSet(), constraint)
                        coefficients = [term.coefficient for term in func.terms]
                        if sum(coefficients .== coefficients[1]) == num_vars # knapsack
                            if s == :lessthan
                                if coefficients[1] > 0
                                    type = QBO_KNAPSACK_UNIT
                                    N = set.upper / coefficients[1]
                                elseif coefficients[1] < 0
                                    type = QBO_KNAPSACK_UNIT_REV
                                    N = set.upper / coefficients[1]
                                end
                            elseif s == :greaterthan
                                if coefficients[1] < 0
                                    type = QBO_KNAPSACK_UNIT
                                    N = set.lower / coefficients[1]
                                elseif coefficients[1] > 0
                                    type = QBO_KNAPSACK_UNIT_REV
                                    N = set.lower / coefficients[1]
                                end
                            elseif s == :equalto
                                type = QBO_KNAPSACK_PROB
                                N = set.value / coefficients[1]
                            end
                        end
                    end
                end
            end
        end

    elseif !isempty(quad_objective_data[1]) # quadratic mixed problem
        type = QMO

    else # linear mixed problem
        type = LMO
    end

    # get the indices of integer variables
    int_vars = get_integer_variables_indices(original_model)

    if heur_params.use_mccormick
        # build the bounds for all models
        global_variable_bounds = build_global_bounds(mccormick_linear_model, heur_params.bigM)
        _ = build_global_bounds(mccormick_original_model, heur_params.bigM)
        _ = build_global_bounds(linear_model, heur_params.bigM)
        _ = build_global_bounds(original_model, heur_params.bigM)
    else
        # build the bounds for all models
        global_variable_bounds = build_global_bounds(linear_model, heur_params.bigM)
        _ = build_global_bounds(original_model, heur_params.bigM)
    end

    if heur_params.use_propagation_global
        model = heur_params.use_mccormick ? mccormick_linear_model : linear_model
        if heur_params.use_external_propagation
            lb_prop, ub_prop, infeasible, num_fixings = apply_prop_external(model, global_variable_bounds, [], fill(0.0, num_variables(model)))
        else
            lb_prop, ub_prop, infeasible, num_fixings = propagate_all(model, global_variable_bounds.lower_bounds, global_variable_bounds.upper_bounds, num_variables(model))
        end
        # TODO We might be infeasible if the bounds we set for continous variables are too tight. In this case we should relax them.
        if !infeasible && num_fixings > 0
            println("Global propagation reduces the domain of $num_fixings variables")
            if heur_params.use_mccormick
                # set bounds for both the linear and the original model
                set_new_bounds!(mccormick_linear_model, lb_prop, ub_prop)
                set_new_bounds!(mccormick_original_model, lb_prop, ub_prop)
            end
            set_new_bounds!(linear_model, lb_prop, ub_prop)
            set_new_bounds!(original_model, lb_prop, ub_prop)
            # update the global bounds
            global_variable_bounds.lower_bounds = lb_prop
            global_variable_bounds.upper_bounds = ub_prop
        end
    end

    cons_lb_list = get_lower_bound_list(original_model)
    cons_ub_list = get_upper_bound_list(original_model)

    @assert num_vars == length(cons_lb_list) == length(cons_ub_list)

    interaction_matrix = create_sparse_quadratic_interaction_matrix(num_variables(original_model), quad_objective_data, quadratic_lowerthan, quadratic_equalto)

    println("\nInstance name: $name")

    write_result_header(result_path, walltime_start)

    # compute diverse vertex covers instead of just one minimum cover
    covers = diverse_vertex_covers(interaction_matrix, heur_params.solver_choice;
                                  num_covers=heur_params.num_covers,
                                  max_overlap_ratio=0.7, time_limit=5.0)

    # the first cover is the min_cover
    min_cover = isempty(covers) ? Int[] : covers[1]

    # if the min cover is not empty, we will use the undercover heuristic
    heur_params.prob_undercover = sizeof(min_cover) > 0 ? heur_params.prob_undercover : 0.0

    if heur_params.mu_barrier == 0.0
        heur_params.mu_barrier = mu_estimate(quad_objective_data)
    end
    @info "Estimated mu: $(heur_params.mu_barrier)"

    initialise_strategy!(heur_params, type, heur_params.mu_barrier)

    if !isempty(result_path) && Threads.nthreads() > 1
        ready_signal = Threads.Event()  # Signal for when the watcher is ready
        # Shared flag for stopping the watcher
        stop_signal = Threads.Atomic{Bool}(false)
        # Start file watcher in a separate thread
        watcher_task = Threads.@spawn file_watcher(result_path, sense == MOI.MIN_SENSE, stop_signal, ready_signal)
        wait(ready_signal) # Ensure the watcher is running before proceeding
    end

    res = nothing

    # Get the number of threads we want to run the heuristic
    if heur_params.num_threads == 0
        # use all available threads
        nthreads = nthreads_for_type()
    else
        # use at most the number of threads available
        nthreads = min(heur_params.num_threads, nthreads_for_type())
    end

    available_threads = nthreads_for_type()

    @info "Threads available: $(available_threads)"
    @info "Starting heuristic with $nthreads threads"
    min_val = Threads.Atomic{Float64}(Inf)

    # First allocation strategy - divide evenly
    threads_per_solver = div(available_threads, nthreads)

    # If we have leftover threads, redistribute them
    remaining_threads = available_threads - (threads_per_solver * nthreads)

    # Create a distribution array where each worker gets its thread allocation
    thread_distribution = fill(threads_per_solver, nthreads)

    # Distribute remaining threads (if any) to maximize usage
    for i in 1:remaining_threads
        thread_distribution[i] += 1
    end

    @info "Thread distribution across $nthreads workers: $thread_distribution (total: $nthreads)"

    Threads.@threads for i in 1:nthreads

        thread_id = Threads.threadid() # WARNING, this is not necessarily i

        # Get thread allocation for this worker
        worker_threads = thread_distribution[i]

        result_file = "ResultFile_t$thread_id"
        heur_stats = HeurStats()
        sol_data = SolutionData(get_best_known_solution("$name"), initial_obj_value, Float64[], 0, result_path, result_file)

        # copies heur_params and changes relevant fields according to the type-specific parallelisation strategy
        heur_params_thread = copy(heur_params)

        parallelisation_strategy!(heur_params_thread, type, spectrum, i, nthreads)

        # Create the optimizer (locally to avoid interference between threads)
        if heur_params.solver_choice == "Gurobi"
            optimizer_original = Gurobi.Optimizer()
            optimizer_linear = Gurobi.Optimizer()
        else
            error("Unsupported solver choice: $(heur_params.solver_choice)")
        end

        # Set number of threads for each solver based on worker's allocation
        MOI.set(optimizer_original, MOI.NumberOfThreads(), worker_threads)
        MOI.set(optimizer_linear, MOI.NumberOfThreads(), worker_threads)
        MOI.set(optimizer_original, MOI.Silent(), !heur_params.verbose_mip_lns)
        MOI.set(optimizer_original, MOI.TimeLimitSec(), heur_params.max_time_lmo)
        MOI.set(optimizer_original, MOI.RelativeGapTolerance(), 2e-4)
        MOI.set(optimizer_linear, MOI.Silent(), !heur_params.verbose_mip_oracle)
        MOI.set(optimizer_linear, MOI.TimeLimitSec(), heur_params.max_time_lmo)

        if heur_params_thread.use_mccormick
            MOI.copy_to(optimizer_original, mccormick_original_model)
            MOI.copy_to(optimizer_linear, mccormick_linear_model)
            problem_data = ProblemData{Tuple{SparseArrays.SparseMatrixCSC{Float64,Int64},Vector{Float64},Float64}}(name,
                optimizer_original, quadratic_lowerthan_mccormick, quadratic_equalto_mccormick, quad_objective_data_mccormick,
                num_added_mccormick_vars, global_variable_bounds, int_vars, cons_lb_list, cons_ub_list, sense, sol_data, time_ref, min_cover,
                covers, (ind1, ind2), N, spectrum, heur_stats, x -> postsolve_mccormick(postsolve_function(x)), type)
        else
            MOI.copy_to(optimizer_original, original_model)
            MOI.copy_to(optimizer_linear, linear_model)
            problem_data = ProblemData{Tuple{SparseArrays.SparseMatrixCSC{Float64,Int64},Vector{Float64},Float64}}(name,
                optimizer_original, quadratic_lowerthan, quadratic_equalto, quad_objective_data,
                0, global_variable_bounds, int_vars, cons_lb_list, cons_ub_list, sense, sol_data, time_ref, min_cover,
                covers, (ind1, ind2), N, spectrum, heur_stats, postsolve_function, type)
        end

        # also this may be good to separate heur_stats from the rest in a cleaner way
        problem_statistics(problem_data, heur_params, i == 1 && heur_params.verbose_statistics)

        try
            # Run the heuristic
            run_heuristic(problem_data, heur_params_thread, optimizer_linear, i, nthreads)
        catch e
            @info "Error in thread $i: $e"
        end

        if i == 1 && heur_params.verbose_statistics
            print_stats(problem_data)
        end

        val = problem_data.sense == MOI.MAX_SENSE ? -problem_data.sol_data.best_found_val : problem_data.sol_data.best_found_val
        Threads.atomic_min!(min_val, val)
    end

    if !isempty(result_path) && Threads.nthreads() > 1
        # Signal the watcher to stop
        stop_signal[] = true
        # Ensure the watcher thread exits cleanly
        wait(watcher_task)
        println("Main computation finished.")
    end

    return min_val[]
end

end # module

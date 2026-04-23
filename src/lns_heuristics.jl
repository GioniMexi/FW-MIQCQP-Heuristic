function copy_global_bounds(global_variable_bounds)
    return VariableBounds(copy(global_variable_bounds.lower_bounds), copy(global_variable_bounds.upper_bounds))
end

function update_heuristic_stats(problem_data, heur_name, obj_val)
    # Determine the objective sense
    objective_sense = MOI.get(problem_data.optimizer, MOI.ObjectiveSense())
    is_better = false

    if objective_sense == MOI.MIN_SENSE
        is_better = obj_val < problem_data.sol_data.best_found_val - 1e-9
    elseif objective_sense == MOI.MAX_SENSE
        is_better = obj_val > problem_data.sol_data.best_found_val + 1e-9
    else
        error("Unknown objective sense")
    end

    if heur_name == "RENS"
        problem_data.heur_stats.num_sols_rens += 1
        if is_better
            problem_data.heur_stats.num_sols_improved_rens += 1
        end
    elseif heur_name == "ARENS"
        problem_data.heur_stats.num_sols_arens += 1
        if is_better
            problem_data.heur_stats.num_sols_improved_arens += 1
        end
    elseif heur_name == "RINS"
        problem_data.heur_stats.num_sols_rins += 1
        if is_better
            problem_data.heur_stats.num_sols_improved_rins += 1
        end
    elseif heur_name == "Undercover"
        problem_data.heur_stats.num_sols_undercover += 1
        if is_better
            problem_data.heur_stats.num_sols_improved_undercover += 1
        end
    elseif heur_name == "Alternating"
        problem_data.heur_stats.num_sols_alternating += 1
        if is_better
            problem_data.heur_stats.num_sols_improved_alternating += 1
        end
    end
end

function optimize_and_get_solution(problem_data, time_limit, heur_name, verbose_statistics)
    time_limit = min(time_limit, 10.0)
    MOI.set(problem_data.optimizer, MOI.TimeLimitSec(), time_limit)
    MOI.optimize!(problem_data.optimizer)
    status = MOI.get(problem_data.optimizer, MOI.TerminationStatus())

    nsols = MOI.get(problem_data.optimizer, MOI.ResultCount())
    if nsols > 0
        xv = [MOI.get(problem_data.optimizer, MOI.VariablePrimal(1), xi) for xi in MOI.get(problem_data.optimizer, MOI.ListOfVariableIndices())]
        obj_val = MOI.get(problem_data.optimizer, MOI.ObjectiveValue())
        if verbose_statistics
            update_heuristic_stats(problem_data, heur_name, obj_val)
        end
        return Vector{Float64}[xv], false
    end

    return Vector{Float64}[], false
end

function rens_bounds(x, problem_data, new_bounds, tol = 1e-5)
    number_fixed_var = 0
    for j in problem_data.int_vars
        new_bounds.lower_bounds[j] = max(new_bounds.lower_bounds[j], floor(x[j]))
        new_bounds.upper_bounds[j] = min(new_bounds.upper_bounds[j], ceil(x[j]))
        if new_bounds.lower_bounds[j] == new_bounds.upper_bounds[j]
            number_fixed_var += 1
        end
        @assert new_bounds.lower_bounds[j] >= problem_data.global_variable_bounds.lower_bounds[j] - tol
        @assert new_bounds.upper_bounds[j] <= problem_data.global_variable_bounds.upper_bounds[j] + tol
    end
    number_of_updated_bounds = update_variable_bounds(problem_data, new_bounds.lower_bounds, new_bounds.upper_bounds)
    return number_fixed_var, number_of_updated_bounds
end

"""
RENS Heuristic:
This heuristic fixes integer variable bounds based solely on the current solution x
by updating the global variable bounds. It then reoptimizes the model within the
remaining time limit.
"""
function rens_heuristic(tree, lmo, x, problem_data, heur_params)
    if heur_params.verbose_lns
        println("\nRunning RENS:")
    end

    new_bounds = copy_global_bounds(problem_data.global_variable_bounds)

    number_fixed_var, number_of_updated_bounds = rens_bounds(x, problem_data, new_bounds, heur_params.eps_int)

    time_limit = heur_params.time_limit - elapsed_time(problem_data.time_ref)

    num_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars

    fixed_int_perc = round(number_fixed_var / num_vars; digits = 2)
    updated_int_bounds_perc = round(number_of_updated_bounds / (2 * num_vars); digits = 2)

    if heur_params.verbose_lns
        println("  RENS percentage fixed: ", fixed_int_perc)
        println("  RENS percentage updated bounds: ", updated_int_bounds_perc)
    end
    # TODO add parameters
    if fixed_int_perc >= heur_params.rens_rins_perc_fixed_vars && time_limit > 2.0
        return optimize_and_get_solution(problem_data, time_limit - 1.0, "RENS", heur_params.verbose_statistics)
    end
    return Vector{Float64}[], false
end

function rins_bounds(x, problem_data, new_bounds, tol = 1e-5)
    number_fixed_var = 0
    for j in problem_data.int_vars
        if x[j] == problem_data.sol_data.best_found_sol[j]
            new_bounds.lower_bounds[j] = problem_data.sol_data.best_found_sol[j]
            new_bounds.upper_bounds[j] = problem_data.sol_data.best_found_sol[j]
            number_fixed_var += 1
        else
            new_bounds.lower_bounds[j] = problem_data.global_variable_bounds.lower_bounds[j]
            new_bounds.upper_bounds[j] = problem_data.global_variable_bounds.upper_bounds[j]
        end
        @assert new_bounds.lower_bounds[j] >= problem_data.global_variable_bounds.lower_bounds[j] - tol
        @assert new_bounds.upper_bounds[j] <= problem_data.global_variable_bounds.upper_bounds[j] + tol
    end
    number_of_updated_bounds = update_variable_bounds(problem_data, new_bounds.lower_bounds, new_bounds.upper_bounds)
    return number_fixed_var, number_of_updated_bounds
end

"""
RINS Heuristic:
This heuristic fixes integer bounds for variables that match the best found solution,
based on the current solution x, and then reoptimizes within the remaining time limit.
"""
function rins_heuristic(tree, lmo, x, problem_data, heur_params)
    if problem_data.sol_data.num_sols > 0
        if heur_params.verbose_lns
            println("\nRunning RINS:")
        end

        new_bounds = copy_global_bounds(problem_data.global_variable_bounds)

        number_fixed_var, number_of_updated_bounds = rins_bounds(x, problem_data, new_bounds, heur_params.eps_int)

        time_limit = heur_params.time_limit - elapsed_time(problem_data.time_ref)

        num_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars

        fixed_int_perc = round(number_fixed_var / num_vars; digits = 2)
        updated_int_bounds_perc = round(number_of_updated_bounds / (2 * num_vars); digits = 2)

        if heur_params.verbose_lns
            println("  RINS percentage fixed: ", fixed_int_perc)
            println("  RINS percentage updated bounds: ", updated_int_bounds_perc)
        end

        # TODO add parameters
        if fixed_int_perc >= heur_params.rens_rins_perc_fixed_vars && time_limit > 2.0
            return optimize_and_get_solution(problem_data, time_limit - 1.0, "RINS", heur_params.verbose_statistics)
        end
    end
    return Vector{Float64}[], false
end

function undercover_prop_bounds(x, problem_data, heur_params, new_bounds, cover_indices=nothing)
    indices = cover_indices === nothing ? problem_data.cover_indices : cover_indices

    if heur_params.use_external_propagation
        new_bounds.lower_bounds, new_bounds.upper_bounds, infeasible, num_reduced_domains = apply_prop_external(problem_data.optimizer, problem_data.global_variable_bounds, indices, x)
    else
        new_bounds.lower_bounds, new_bounds.upper_bounds, infeasible, num_reduced_domains = apply_prop(problem_data, indices, x)
    end
    if infeasible
        return 0, 0, infeasible, 0
    end
    number_of_updated_bounds = update_variable_bounds(problem_data, new_bounds.lower_bounds, new_bounds.upper_bounds)
    return length(indices), number_of_updated_bounds, false, num_reduced_domains
end

function undercover_bounds(x, problem_data, new_bounds, cover_indices=nothing)
    tol = 1e-5
    indices = cover_indices === nothing ? problem_data.cover_indices : cover_indices

    number_fixed_var = 0

    for j in indices
        if j in problem_data.int_vars
            new_bounds.lower_bounds[j] = round(x[j])
        else
            new_bounds.lower_bounds[j] = x[j]
        end
        new_bounds.upper_bounds[j] = new_bounds.lower_bounds[j]
        number_fixed_var += 1
        @assert new_bounds.lower_bounds[j] >= problem_data.global_variable_bounds.lower_bounds[j] - tol
        @assert new_bounds.upper_bounds[j] <= problem_data.global_variable_bounds.upper_bounds[j] + tol
    end

    number_of_updated_bounds = update_variable_bounds(problem_data, new_bounds.lower_bounds, new_bounds.upper_bounds)
    return number_fixed_var, number_of_updated_bounds, false
end

"""
Undercover Heuristic:
This heuristic updates variable bounds by propagating coverings computed from the
current solution x. It then reoptimizes the model within the remaining time limit.
Iterates through multiple covers if available.
"""
function undercover_heuristic(tree, lmo, x, problem_data, heur_params)
    if heur_params.verbose_lns
        println("\nRunning Undercover:")
    end

    covers_to_try = isempty(problem_data.all_covers) ?
                    [problem_data.cover_indices] :
                    problem_data.all_covers

    for (cover_idx, current_cover) in enumerate(covers_to_try)
        if heur_params.verbose_lns
            println("  Trying cover #$cover_idx of $(length(covers_to_try))")
        end

        new_bounds = copy_global_bounds(problem_data.global_variable_bounds)
        if heur_params.use_propagation_undercover
            number_fixed_var, number_of_updated_bounds, infeasible, _ = undercover_prop_bounds(x, problem_data, heur_params, new_bounds, current_cover)
        else
            number_fixed_var, number_of_updated_bounds, infeasible = undercover_bounds(x, problem_data, new_bounds, current_cover)
        end

        if infeasible
            if heur_params.verbose_lns
                println("  Cover #$cover_idx led to infeasibility, trying next cover")
            end
            continue
        end

        time_limit = heur_params.time_limit - elapsed_time(problem_data.time_ref)

        num_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars

        fixed_int_perc = round(number_fixed_var / num_vars; digits = 2)
        updated_int_bounds_perc = round(number_of_updated_bounds / (2 * num_vars); digits = 2)

        if heur_params.verbose_lns
            println("  Undercover cover #$cover_idx size: $(length(current_cover))")
            println("  Undercover percentage fixed: ", fixed_int_perc)
            println("  Undercover percentage updated bounds: ", updated_int_bounds_perc)
        end

        # Try to solve with this cover
        if time_limit > 2.0
            sols, early_term = optimize_and_get_solution(problem_data, time_limit - 1.0, "Undercover", heur_params.verbose_statistics)
            if heur_params.verbose_lns
                println("  Undercover cover #$cover_idx led to $(length(sols)) solutions")
            end
            # If we found a solution with this cover, return it
            if !isempty(sols)
                return sols, early_term
            end
        end
    end

    return Vector{Vector{Float64}}(), false
end

function active_set_rens_bounds(x, problem_data, active_set, new_bounds, tol = 1e-5)
    number_original_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars
    number_fixed_var = 0
    for j in get_all_variables_indices(problem_data.optimizer)
        if j > number_original_vars
            continue
        end
        lb_j = new_bounds.upper_bounds[j]
        ub_j = new_bounds.lower_bounds[j]
        for elem in active_set
            lb_j = min(lb_j, elem[2][j])
            ub_j = max(ub_j, elem[2][j])
        end
        # TODO: possibly change domains of continuous variables
        if j in problem_data.int_vars
            new_bounds.lower_bounds[j] = min(floor(x[j] + tol), lb_j)
            new_bounds.upper_bounds[j] = max(ceil(x[j] - tol), ub_j)
        else
            new_bounds.lower_bounds[j] = min(x[j], lb_j)
            new_bounds.upper_bounds[j] = max(x[j], ub_j)
        end

        # make sure the new bounds are within the global bounds
        new_bounds.lower_bounds[j] = max(new_bounds.lower_bounds[j], problem_data.global_variable_bounds.lower_bounds[j])
        new_bounds.upper_bounds[j] = min(new_bounds.upper_bounds[j], problem_data.global_variable_bounds.upper_bounds[j])

        if new_bounds.lower_bounds[j] == new_bounds.upper_bounds[j]
            number_fixed_var += 1
        end
    end
    number_of_updated_bounds = update_variable_bounds(problem_data, new_bounds.lower_bounds, new_bounds.upper_bounds)
    return number_fixed_var, number_of_updated_bounds
end

"""
Active set RENS Heuristic:
This heuristic fixes integer variable bounds based on the active set
and the current solution x.
It then reoptimizes the model within the remaining time limit.
"""
function active_set_rens_heuristic(tree, lmo, x, problem_data, heur_params)
    if heur_params.verbose_lns
        println("\nRunning Active Set RENS:")
    end

    node = tree.nodes[tree.root.current_node_id[]]

    if length(node.active_set) >= 2
        new_bounds = copy_global_bounds(problem_data.global_variable_bounds)

        number_fixed_var, number_of_updated_bounds = active_set_rens_bounds(x, problem_data, node.active_set, new_bounds, heur_params.eps_int)

        time_limit = heur_params.time_limit - elapsed_time(problem_data.time_ref)

        num_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars

        fixed_int_perc = round(number_fixed_var / num_vars; digits = 2)
        updated_int_bounds_perc = round(number_of_updated_bounds / (2 * num_vars); digits = 2)

        if heur_params.verbose_lns
            println("  Active set points used: ", length(node.active_set))
            println("  Active Set RENS percentage fixed: ", fixed_int_perc)
            println("  Active Set RENS percentage updated bounds: ", updated_int_bounds_perc)
        end
        # TODO add parameters
        if fixed_int_perc >= heur_params.rens_rins_perc_fixed_vars && time_limit > 2.0
            return optimize_and_get_solution(problem_data, time_limit - 1.0, "ARENS", heur_params.verbose_statistics)
        end
    end
    return Vector{Float64}[], false
end

"""
Alternating heuristic for all points in the active set. Only for QUBOs with bipartite structure.
"""
function alternating_heuristic(tree, lmo, x, problem_data, heur_params)
    node = tree.nodes[tree.root.current_node_id[]]

    if length(node.active_set) >= 1
        if heur_params.verbose_lns
            println("\nRunning Alternating:")
        end

        Q, b, c0 = problem_data.quad_objective_data
        ind1, ind2 = problem_data.bipartition_indices
        if isempty(ind1) || isempty(ind2)
            return Vector{Float64}[], false
        end

        # x1 and x2 encode the variables via their sign: negative means activated variable
        # they are modified in place in the alternating_min auxiliary functions
        x1 = zeros(length(ind1))
        x2 = zeros(length(ind2))
        best_obj = Inf
        best_x1 = zeros(length(ind1))
        best_x2 = zeros(length(ind2))

        for (w, a) in node.active_set
            for i in eachindex(x1)
                x1[i] = a[ind1[i]] ≈ 1 ? -1 : 0
            end
            obj = alternating_min1!(x1, x2, Q, b, ind1, ind2)
            if obj < best_obj
                best_obj = obj
                best_x1 .= x1
                best_x2 .= x2
            end
            for j in eachindex(x2)
                x2[j] = a[ind2[j]] ≈ 1 ? -1 : 0
            end
            obj = alternating_min2!(x1, x2, Q, b, ind1, ind2)
            if obj < best_obj
                best_obj = obj
                best_x1 .= x1
                best_x2 .= x2
            end
        end

        xv = zeros(length(ind1) + length(ind2))
        @views xv[ind1] .= best_x1 .< 0
        @views xv[ind2] .= best_x2 .< 0
        if heur_params.verbose_statistics
            update_heuristic_stats(problem_data, "Alternating", best_obj + c0)
        end
        return Vector{Float64}[xv], false
    end

    return Vector{Float64}[], false
end

# Alternating maximisation starting from a given x1
function alternating_min1!(x1, x2, Q, b, ind1, ind2)
    obj = 0.0
    obj_last = NaN
    while obj != obj_last
        obj_last = obj
        @views x2 .= b[ind2]
        for i in eachindex(x1)
            if x1[i] < 0 # corresponding variable is activated
                @views x2 .+= Q[ind1[i], ind2]
            end
        end
        obj = 0
        @views x1 .= b[ind1]
        for j in eachindex(x2)
            if x2[j] < 0 # corresponding variable is activated
                obj += b[ind2[j]]
                @views x1 .+= Q[ind1, ind2[j]]
            end
        end
        for i in eachindex(x1)
            if x1[i] < 0
                obj += x1[i]
            end
        end
    end
    return obj
end

# Alternating maximisation starting from a given x2
function alternating_min2!(x1, x2, Q, b, ind1, ind2)
    obj = 0.0
    obj_last = NaN
    while obj != obj_last
        obj_last = obj
        @views x1 .= b[ind1]
        for j in eachindex(x2)
            if x2[j] < 0 # corresponding variable is activated
                @views x1 .+= Q[ind1, ind2[j]]
            end
        end
        obj = 0
        @views x2 .= b[ind2]
        for i in eachindex(x1)
            if x1[i] < 0 # corresponding variable is activated
                obj += b[ind1[i]]
                @views x2 .+= Q[ind1[i], ind2]
            end
        end
        for j in eachindex(x2)
            if x2[j] < 0
                obj += x2[j]
            end
        end
    end
    return obj
end

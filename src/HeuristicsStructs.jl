"""
Container for tracking heuristic statistics.
"""
mutable struct HeurStats
    num_vars::Int64
    num_quad_cons::Int64
    num_bin_vars::Int64
    num_int_vars::Int64
    num_cont_vars::Int64
    num_quad_eq_cons::Int64
    num_quad_leq_cons::Int64
    num_quad_non_bilinear::Int64
    num_quad_bilinear::Int64
    num_quad_non_convex::Int64
    num_quad_convex::Int64
    num_sols_rens::Int64
    num_sols_arens::Int64
    num_sols_rins::Int64
    num_sols_undercover::Int64
    num_sols_alternating::Int64
    num_sols_improved_rens::Int64
    num_sols_improved_arens::Int64
    num_sols_improved_rins::Int64
    num_sols_improved_undercover::Int64
    num_sols_improved_alternating::Int64

    # Constructor with default values set to zero
    HeurStats() = new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

"""
Container for the heuristic parameters.
"""
mutable struct HeurParams
    solver_choice::String
    # Competition
    time_limit::Float64
    node_limit::Int
    eps_cons::Float64
    eps_int::Float64
    # Propagation
    bigM::Float64
    bigM_McCormick::Float64
    # Quadratics Handling
    shift_barrier::Int
    mu_barrier::Float64
    power_penalty::Float64 # -1 to use the log barrier instead
    adaptive_mu::Bool      # Whether to adapt mu dynamically based on solutions
    # LNS Heuristics
    prob_rens::Float64
    prob_arens::Float64
    prob_rins::Float64
    prob_undercover::Float64
    prob_alternating::Float64
    rens_rins_perc_fixed_vars::Float64
    use_mccormick::Bool
    use_propagation_global::Bool
    use_propagation_undercover::Bool
    use_external_propagation::Bool
    num_covers::Int
    # Custom LMO
    prob_custom_heu::Float64
    # Boscia heuristics
    prob_follow_gradient_heu::Float64
    max_iter_follow_gradient_heu::Int
    prob_rounding_01_heu::Float64
    prob_probability_rounding_heu::Float64
    # Boscia settings
    use_secant::Bool
    fw_lazy::Bool
    fw_max_iter::Int
    fw_min_iter::Int
    fw_variant::Int
    fw_timeout::Float64
    # Verbosity
    verbose_statistics::Bool
    verbose_problem_structure_statistics::Bool
    verbose_boscia::Bool
    verbose_fw::Bool
    verbose_mip_oracle::Bool
    verbose_mip_lns::Bool
    verbose_violation::Bool
    verbose_lns::Bool
    # Convexification
    convexify_objective::Float64
    deperspective::Bool
    # Restart
    restart::Bool
    restart_min_time::Int
    # fine tuning
    perform_benchmark_oracle::Bool
    num_threads::Int
    max_time_lmo::Float64
end

function HeurParams(;
    solver_choice = "Gurobi",
    time_limit = 300.0,
    node_limit = 100,
    eps_cons = 1e-6,
    eps_int = 1e-5,
    bigM = 1000000.0,
    bigM_McCormick = 10000.0,
    shift_barrier = 1,
    mu_barrier = 0.0,
    power_penalty = 1.2,
    adaptive_mu = false,
    prob_rens = 0.0,
    prob_arens = 1.0,
    prob_rins = 1.0,
    prob_undercover = 1.0,
    prob_alternating = 0.0,
    rens_rins_perc_fixed_vars = 0.5,
    use_mccormick = false,
    use_propagation_global = true,
    use_propagation_undercover = false,
    use_external_propagation = true,
    num_covers = 1,
    prob_custom_heu = 0.0,
    prob_follow_gradient_heu = 1.0,
    max_iter_follow_gradient_heu = 100,
    prob_rounding_01_heu = 1.0,
    prob_probability_rounding_heu = 1.0,
    use_secant = true,
    fw_lazy = true,
    fw_max_iter = 500,
    fw_min_iter = 50,
    fw_variant = 0, # 0: default, 1: vanilla
    fw_timeout = 15.0,
    verbose_statistics = false,
    verbose_problem_structure_statistics = false,
    verbose_boscia = false,
    verbose_fw = false,
    verbose_mip_oracle = false,
    verbose_mip_lns = false,
    verbose_violation = false,
    verbose_lns = false,
    convexify_objective = 0.0,
    deperspective = true,
    restart = true,
    restart_min_time = 10,
    perform_benchmark_oracle = false,
    num_threads = 7,
    max_time_lmo = 10.0
)
    return HeurParams(
        solver_choice, time_limit, node_limit, eps_cons, eps_int, bigM, bigM_McCormick, shift_barrier, mu_barrier,
        power_penalty, adaptive_mu, prob_rens, prob_arens, prob_rins, prob_undercover, prob_alternating, rens_rins_perc_fixed_vars,
        use_mccormick, use_propagation_global, use_propagation_undercover, use_external_propagation, num_covers, prob_custom_heu,
        prob_follow_gradient_heu, max_iter_follow_gradient_heu, prob_rounding_01_heu, prob_probability_rounding_heu, use_secant, fw_lazy, fw_max_iter, fw_min_iter, fw_variant,
        fw_timeout, verbose_statistics, verbose_problem_structure_statistics, verbose_boscia, verbose_fw, verbose_mip_oracle, verbose_mip_lns, verbose_violation,
        verbose_lns, convexify_objective, deperspective, restart, restart_min_time, perform_benchmark_oracle, num_threads, max_time_lmo
    )
end

function Base.copy(h::HeurParams)
    return HeurParams(
        h.solver_choice, h.time_limit, h.node_limit, h.eps_cons, h.eps_int, h.bigM, h.bigM_McCormick, h.shift_barrier, h.mu_barrier,
        h.power_penalty, h.adaptive_mu, h.prob_rens, h.prob_arens, h.prob_rins, h.prob_undercover, h.prob_alternating, h.rens_rins_perc_fixed_vars,
        h.use_mccormick, h.use_propagation_global, h.use_propagation_undercover, h.use_external_propagation, h.num_covers, h.prob_custom_heu,
        h.prob_follow_gradient_heu, h.max_iter_follow_gradient_heu, h.prob_rounding_01_heu, h.prob_probability_rounding_heu, h.use_secant, h.fw_lazy, h.fw_max_iter ,h.fw_min_iter, h.fw_variant,
        h.fw_timeout, h.verbose_statistics, h.verbose_problem_structure_statistics, h.verbose_boscia, h.verbose_fw, h.verbose_mip_oracle, h.verbose_mip_lns, h.verbose_violation,
        h.verbose_lns, h.convexify_objective, h.deperspective, h.restart, h.restart_min_time, h.perform_benchmark_oracle, h.num_threads, h.max_time_lmo
    )
end

"""
Container for the solution data.
"""
mutable struct SolutionData
    best_known_val::Float64
    best_found_val::Float64
    best_found_sol::Vector{Float64}
    num_sols::Int64
    result_path::String
    result_file::String
end

"""
Enum for the type of problem.
"""
@enum QType::Int8 begin
    UNKNOWN = 0
    LBO = 10                   # LBO: linear binary optimization
    QUBO = 11                  # QUBO: quadratic unconstrained binary optimization
    QUBO_BIPARTITE = 12        #       bipartite underlying graph for alternating heuristic
    QBO = 13                   # QBO: quadratic binary optimization
    QBO_KNAPSACK_PROB = 14     #      with knapsack equality constraint
    QBO_KNAPSACK_UNIT = 15     #      with knapsack inequality constraint
    QBO_KNAPSACK_UNIT_REV = 16 #      with reversed knapsack inequality constraint
    QMO = 20                   # QMO: quadratic mixed optimization
    LMO = 30                   # LMO: linear mixed optimization
end

function is_type_lmo(type::QType)
    return type in [LBO, LMO]
end

function is_type_qbo(type::QType)
    return type in [QUBO, QUBO_BIPARTITE, QBO, QBO_KNAPSACK_PROB, QBO_KNAPSACK_UNIT, QBO_KNAPSACK_UNIT_REV]
end

struct ProblemData{T}
    prob_name::String
    optimizer::MOI.AbstractOptimizer
    quadratic_lowerthan::Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}},Tuple}
    quadratic_equalto::Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}},Tuple}
    quad_objective_data::T
    num_auxiliary_vars::Int64
    global_variable_bounds::VariableBounds
    int_vars::Vector{Int64}
    cons_lb_list::Vector{MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}}
    cons_ub_list::Vector{MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}}
    sense::MOI.OptimizationSense
    sol_data::SolutionData
    time_ref::DateTime
    cover_indices::Vector{Int64}
    all_covers::Vector{Vector{Int64}}
    bipartition_indices::NTuple{2,Vector{Int64}}
    N::Float64
    spectrum::Vector{Float64}
    heur_stats::HeurStats
    postsolve_function::Function
    type::QType
end

# Define a method for the show function to print all fields of HeurParams
function Base.show(io::IO, params::HeurParams)
    println(io, "Parameters:")
    for field in fieldnames(HeurParams)
        println(io, "  $(field): $(getfield(params, field))")
    end
end

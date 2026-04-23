"""
    TimeTrackingLMO  <: FW.LMO

An LMO wrapping another one tracking the time, number of nodes and number of calls.

"""
mutable struct TimeTrackingLMO{BLMO<:BoundedLinearMinimizationOracle, D<:Dates.DateTime} <:
               FrankWolfe.LinearMinimizationOracle
    blmo::BLMO
    optimizing_times::Vector{Float64}
    optimizing_nodes::Vector{Int}
    simplex_iterations::Vector{Int}
    ncalls::Int
    int_vars::Vector{Int}
    time_ref::D
    type_moi::Bool
    time_limit::Float64
    max_time_lmo::Float64
end

TimeTrackingLMO(blmo::BoundedLinearMinimizationOracle, time_ref, time_limit, max_time_lmo) =
    TimeTrackingLMO(blmo, Float64[], Int[], Int[], 0, Int[], time_ref, isa(blmo, MathOptBLMO), time_limit, max_time_lmo)

TimeTrackingLMO(blmo::BoundedLinearMinimizationOracle, int_vars, time_ref, time_limit, max_time_lmo) =
    TimeTrackingLMO(blmo, Float64[], Int[], Int[], 0, int_vars, time_ref, isa(blmo, MathOptBLMO), time_limit, max_time_lmo)

# if we want to reset the info between nodes in Bonobo
function reset!(tlmo::TimeTrackingLMO)
    empty!(tlmo.optimizing_times)
    empty!(tlmo.optimizing_nodes)
    empty!(tlmo.simplex_iterations)
    return tlmo.ncalls = 0
end

function FrankWolfe.compute_extreme_point(tlmo::TimeTrackingLMO, d; kwargs...)
    tlmo.ncalls += 1
    free_model(tlmo.blmo)
    if tlmo.type_moi && isfinite(tlmo.time_limit)
        time_limit = tlmo.time_limit - float(Dates.value(Dates.now() - tlmo.time_ref))/1000
        time_limit = time_limit <= 0.1 ? tlmo.max_time_lmo : time_limit
        MOI.set(tlmo.blmo.o, MOI.TimeLimitSec(), min(time_limit, tlmo.max_time_lmo))
    end
    v = FrankWolfe.compute_extreme_point(tlmo.blmo, d; kwargs)

    if !is_linear_feasible(tlmo, v)
        @debug "Vertex not linear feasible $(v)"
        @assert is_linear_feasible(tlmo, v)
    end
    v[tlmo.int_vars] = round.(v[tlmo.int_vars])

    opt_times, numberofnodes, simplex_iterations = get_BLMO_solve_data(tlmo.blmo)

    push!(tlmo.optimizing_times, opt_times)
    push!(tlmo.optimizing_nodes, numberofnodes)
    push!(tlmo.simplex_iterations, simplex_iterations)

    free_model(tlmo.blmo)
    return v
end

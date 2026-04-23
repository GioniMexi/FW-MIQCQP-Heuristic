"""
    VariableBounds

Keeps track of the bounds of the integer (binary) variables.

`lower_bounds` dictionary of Float64, index is the key.
`upper_bounds` dictionary of Float64, index is the key.
"""
mutable struct VariableBounds
    lower_bounds::Dict{Int,Float64}
    upper_bounds::Dict{Int,Float64}
end

VariableBounds() = VariableBounds(Dict{Int,Float64}(), Dict{Int,Float64}())

function Base.push!(ib::VariableBounds, (idx, bound), sense::Symbol)
    if sense == :greaterthan
        ib.lower_bounds[idx] = bound
    elseif sense == :lessthan
        ib.upper_bounds[idx] = bound
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
    return ib
end

function Base.isempty(ib::VariableBounds)
    return isempty(ib.lower_bounds) && isempty(ib.upper_bounds)
end

Base.copy(ib::VariableBounds) = VariableBounds(copy(ib.lower_bounds), copy(ib.upper_bounds))

# convenient call
# ib[3, :lessthan] or ib[3, :greaterthan]
function Base.getindex(ib::VariableBounds, idx::Int, sense::Symbol)
    if sense == :lessthan
        ib.upper_bounds[idx]
    elseif sense == :greaterthan
        ib.lower_bounds[idx]
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

function Base.get(ib::VariableBounds, (idx, sense), default)
    if sense == :lessthan
        get(ib.upper_bounds, idx, default)
    elseif sense == :greaterthan
        get(ib.lower_bounds, idx, default)
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

function Base.setindex!(ib::VariableBounds, val, idx::Int, sense::Symbol)
    if sense == :lessthan
        ib.upper_bounds[idx] = val
    elseif sense == :greaterthan
        ib.lower_bounds[idx] = val
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

function Base.haskey(ib::VariableBounds, (idx, sense))
    if sense == :lessthan
        return haskey(ib.upper_bounds, idx)
    elseif sense == :greaterthan
        return haskey(ib.lower_bounds, idx)
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

"""
Get the list of lower bounds.
"""
function get_lower_bound_list(o)
    return MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{Float64}}())
end

"""
Get the list of upper bounds.
"""
function get_upper_bound_list(o)
    return MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{Float64}}())
end

"""
Read bound value for c_idx.
"""
function get_bound(o, c_idx, sense::Symbol)
    return MOI.get(o, MOI.ConstraintSet(), c_idx)
end

"""
Check if the subject of the bound c_idx is an integer variable (recorded in int_vars).
"""
function is_constraint_on_int_var(c_idx, int_vars)
    return c_idx.value in int_vars
end

"""
To check if there is bound for the variable in the global or node bounds.
"""
function is_bound_in(o, c_idx, bounds)
    return haskey(bounds, c_idx.value)
end

"""
Delete bounds.
"""
function delete_bounds!(o, cons_delete)
    for (d_idx, _) in cons_delete
        MOI.delete(o, d_idx)
    end
end

"""
Add bound constraint.
"""
function add_bound_constraint!(o, key, value, sense::Symbol)
    if sense == :lessthan
        MOI.add_constraint(o, MOI.VariableIndex(key), MOI.LessThan(value))
    elseif sense == :greaterthan
        MOI.add_constraint(o, MOI.VariableIndex(key), MOI.GreaterThan(value))
    end
end

"""
Get the index of the variable the bound is working on.
"""
function get_index_var(c_idx)
    return c_idx.value
end

"""
Change the value of the bound c_idx.
"""
function set_bound!(o, c_idx, value, sense::Symbol)
    if sense == :lessthan
        MOI.set(o, MOI.ConstraintSet(), c_idx, MOI.LessThan(value))
    elseif sense == :greaterthan
        MOI.set(o, MOI.ConstraintSet(), c_idx, MOI.GreaterThan(value))
    else
        error("Allowed values for sense are :lessthan and :greaterthan!")
    end
end

"""
Add explicit bounds for binary variables.
"""
function explicit_bounds_binary_var(o)
    # adding binary bounds explicitly
    binary_variables = get_binary_variables(o)
    for idx in binary_variables
        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx.value)
        if !MOI.is_valid(o, cidx)
            MOI.add_constraint(o, MOI.VariableIndex(idx.value), MOI.LessThan(1.0))
        end
        @assert MOI.is_valid(o, cidx)
        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx.value)
        if !MOI.is_valid(o, cidx)
            MOI.add_constraint(o, MOI.VariableIndex(idx.value), MOI.GreaterThan(0.0))
        end
    end
end

"""
Bad models may contain some fixed variables. This function replaces fixed variables with two bounds.
"""
function replace_fixed_vars!(o, variable_indices)
    for idx in variable_indices
        cidx_eq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}}(idx)
        if MOI.is_valid(o, cidx_eq)
            s = MOI.get(o, MOI.ConstraintSet(), cidx_eq)
            MOI.delete(o, cidx_eq)
            MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.LessThan(s.value))
            MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.GreaterThan(s.value))
        end
    end
end

"""
Bound all variables for the LMO
"""
function bound_variables!(o, variable_indices, bigM)

    # bad models may contain fixed variables
    replace_fixed_vars!(o, variable_indices)

    # add binary bounds explicitly
    explicit_bounds_binary_var(o)

    #TODO we set the bounds to -bigM and bigM if they are not given...
    for idx in variable_indices
        cidx_leq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx)
        cidx_geq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx)
        if !MOI.is_valid(o, cidx_leq) && !MOI.is_valid(o, cidx_geq)
            MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.LessThan(bigM))
            MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.GreaterThan(-bigM))
            # if at least one of the bounds is given
        else
            for ST in (MOI.LessThan{Float64}, MOI.GreaterThan{Float64})
                cidx = MOI.ConstraintIndex{MOI.VariableIndex,ST}(idx)
                if !MOI.is_valid(o, cidx)
                    # if upper bound is not given set it to lower bound + bigM
                    if ST == MOI.LessThan{Float64}
                        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx)
                        @assert MOI.is_valid(o, cidx)
                        s = MOI.get(o, MOI.ConstraintSet(), cidx)
                        MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.LessThan(s.lower + bigM))
                        # if lower bound is not given set it to upper bound - bigM
                    elseif ST == MOI.GreaterThan{Float64}
                        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx)
                        @assert MOI.is_valid(o, cidx)
                        s = MOI.get(o, MOI.ConstraintSet(), cidx)
                        MOI.add_constraint(o, MOI.VariableIndex(idx), MOI.GreaterThan(s.upper - bigM))
                    end
                end
            end
        end
    end
end

"""
Get the bounds for all binary and integer variables.
"""
function build_global_bounds(o, bigM)
    variable_indices = get_all_variables_indices(o)

    global_bounds = VariableBounds()

    # add variable bounds explicity for binary or unbounded variables
    bound_variables!(o, variable_indices, bigM)

    # get the bounds for all variables
    for idx in variable_indices
        cidx_leq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx)
        cidx_geq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx)
        @assert MOI.is_valid(o, cidx_leq) && MOI.is_valid(o, cidx_geq)
        s_leq = MOI.get(o, MOI.ConstraintSet(), cidx_leq)
        s_geq = MOI.get(o, MOI.ConstraintSet(), cidx_geq)
        global_bounds[idx, :lessthan] = s_leq.upper
        global_bounds[idx, :greaterthan] = s_geq.lower
    end
    return global_bounds
end

"""
Update bounds and returns the number of updated bounds.
"""
function update_variable_bounds(problem_data, reference_lb, reference_ub)
    number_original_vars = num_variables(problem_data.optimizer) - problem_data.num_auxiliary_vars
    number_of_updated_bounds = 0
    # Lower bounds on variables
    for c_idx in problem_data.cons_lb_list
        @assert is_bound_in(problem_data.optimizer, c_idx, problem_data.global_variable_bounds.lower_bounds)
        v_idx = get_index_var(c_idx)
        if v_idx > number_original_vars
            continue
        end
        @assert reference_lb[v_idx] >= problem_data.global_variable_bounds.lower_bounds[v_idx] - 1e-5
        @assert reference_lb[v_idx] <= problem_data.global_variable_bounds.upper_bounds[v_idx] + 1e-5
        if reference_lb[v_idx] > problem_data.global_variable_bounds.lower_bounds[v_idx]
            number_of_updated_bounds += 1
        end
        set_bound!(problem_data.optimizer, c_idx, reference_lb[v_idx], :greaterthan)
    end
    # Upper bounds on variables
    for c_idx in problem_data.cons_ub_list
        @assert is_bound_in(problem_data.optimizer, c_idx, problem_data.global_variable_bounds.upper_bounds)
        v_idx = get_index_var(c_idx)
        if v_idx > number_original_vars
            continue
        end
        @assert reference_ub[v_idx] >= problem_data.global_variable_bounds.lower_bounds[v_idx] - 1e-5
        @assert reference_ub[v_idx] <= problem_data.global_variable_bounds.upper_bounds[v_idx] + 1e-5
        if reference_ub[v_idx] < problem_data.global_variable_bounds.upper_bounds[v_idx]
            number_of_updated_bounds += 1
        end
        set_bound!(problem_data.optimizer, c_idx, reference_ub[v_idx], :lessthan)
    end
    return number_of_updated_bounds
end

"""
Set new bounds for the variables.
"""
# TODO: check this!
function set_new_bounds!(optimizer, lower_bounds, upper_bounds)
    num_vars = num_variables(optimizer)
    for (idx, lb) in lower_bounds
        # avoid setting bounds for potentially inexistent McCormick variables
        if idx > num_vars
            continue
        end
        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx)
        @assert MOI.is_valid(optimizer, cidx)
        s = MOI.get(optimizer, MOI.ConstraintSet(), cidx)
        if lb > s.lower + 1e-5
            c_idx_leq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx)
            @assert lb <= MOI.get(optimizer, MOI.ConstraintSet(), c_idx_leq).upper
            set_bound!(optimizer, cidx, lb, :greaterthan)
        end
    end
    for (idx, ub) in upper_bounds
        if idx > num_vars
            continue
        end
        cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(idx)
        @assert MOI.is_valid(optimizer, cidx)
        s = MOI.get(optimizer, MOI.ConstraintSet(), cidx)
        if ub < s.upper - 1e-5
            c_idx_geq = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(idx)
            @assert ub >= MOI.get(optimizer, MOI.ConstraintSet(), c_idx_geq).lower
            set_bound!(optimizer, cidx, ub, :lessthan)
        end
    end
end

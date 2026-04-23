# This file is not part of the module, we still leave it here for testing purposes
using SparseArrays
using MathOptInterface
const MOI = MathOptInterface

const DEBUG_PROPAGATION = false  # Set to true to enable debug output

# Consistent tolerance for propagation (matching C++ implementation)
const PROP_TOL = 1e-6
const PROP_TOL_STRICT = 1e-9
const INF_BOUND = 1e20  # Bound beyond which we consider as infinite

"""
    floor_eps(x)

Floor with epsilon tolerance, matching C++ floorEps behavior.
"""
floor_eps(x::Float64) = floor(x + PROP_TOL_STRICT)

"""
    ceil_eps(x)

Ceiling with epsilon tolerance, matching C++ ceilEps behavior.
"""
ceil_eps(x::Float64) = ceil(x - PROP_TOL_STRICT)

"""
    LinearConstraintData

Holds a scalar affine function and its sparse representation along with
additional data extracted from the constraint.

Fields:
- `f`: The original affine function (`MOI.ScalarAffineFunction{Float64}`).
- `b`: Sparse coefficient vector (`SparseVector{Float64,Int}`) built from `f`.
- `c0`: The constant term extracted from `f`.
- `rhs`: The right-hand side value from the constraint's set.
- `sense`: The constraint type (`:EqualTo`, `:LessThan`, or `:GreaterThan`).
"""
struct LinearConstraintData
    f::MOI.ScalarAffineFunction{Float64}
    b::SparseVector{Float64,Int}
    c0::Float64
    rhs::Float64
    sense::Symbol
end

"""
    PropConstraintData

Augments a `LinearConstraintData` with incremental propagation information.

Fields:
- `lin`: A `LinearConstraintData` instance.
- `min_activity`: The total minimal activity of the constraint (finite part).
- `max_activity`: The total maximal activity of the constraint (finite part).
- `min_act_inf_cnt`: Number of variables contributing -Inf to min activity.
- `max_act_inf_cnt`: Number of variables contributing +Inf to max activity.
- `contributions`: A dictionary mapping variable indices (`Int`) to a tuple
  `(contribution_min, contribution_max)`.
"""
mutable struct PropConstraintData
    lin::LinearConstraintData
    min_activity::Float64
    max_activity::Float64
    min_act_inf_cnt::Int
    max_act_inf_cnt::Int
    contributions::Dict{Int,Tuple{Float64,Float64}}
end

"""
    constraints_for_var(var_to_constraints, var)

Return a vector of constraint indices that mention the variable `var` using the
mapping in `var_to_constraints`. If `var` is not present, an empty vector is returned.
"""
function constraints_for_var(var_to_constraints::Dict{Int,Vector{Int}}, var::Int)::Vector{Int}
    return get(var_to_constraints, var, Int[])
end

"""
    build_linear_function_data_sparse(f, nvars)

Given a `MOI.ScalarAffineFunction{Float64}` `f` and the total number of variables `nvars`,
construct a sparse coefficient vector for the function.

# Arguments
- `f::MOI.ScalarAffineFunction{Float64}`: The affine function.
- `nvars::Int`: The total number of variables (used to set the length of the sparse vector).

# Returns
A tuple `(b, c0)` where:
- `b` is a `SparseVector{Float64,Int}` containing the coefficients for each variable.
- `c0` is the constant term of the function.

Assumes that each term in `f.terms` has a field `variable.value` giving the variable index.
"""
function build_linear_function_data_sparse(f::MOI.ScalarAffineFunction{Float64}, nvars::Int)
    c0 = f.constant
    b = spzeros(nvars)
    for term in f.terms
        b[term.variable.value] = term.coefficient
    end
    return b, c0
end

"""
    get_rhs_and_type(cons_set)

Extracts the right-hand side value and constraint type symbol from a constraint set.

# Arguments
- `cons_set`: A constraint set, which should be one of `MOI.LessThan{Float64}`,
  `MOI.GreaterThan{Float64}`, or `MOI.EqualTo{Float64}`.

# Returns
A tuple `(rhs, sense)` where:
- `rhs` is the numeric right-hand side value.
- `sense` is a symbol: `:LessThan`, `:GreaterThan`, or `:EqualTo`.

If `cons_set` is not one of the expected types, the function throws an error.
"""
function get_rhs_and_type(cons_set)
    if cons_set isa MOI.LessThan{Float64}
        return cons_set.upper, :LessThan
    elseif cons_set isa MOI.GreaterThan{Float64}
        return cons_set.lower, :GreaterThan
    elseif cons_set isa MOI.EqualTo{Float64}
        return cons_set.value, :EqualTo
    else
        error("Unknown constraint set: $(cons_set)")
    end
end

"""
    get_all_linear_constraints(o, nvars)

Extracts all linear constraints from the optimizer `o` and returns them as a vector of
`LinearConstraintData` objects. This function processes equality, less-than, and greater-than
constraints separately.

# Arguments
- `o`: The optimizer object from which constraints are retrieved.
- `nvars::Int`: The total number of variables in the model (used for building the sparse vector).

# Returns
A `Vector{LinearConstraintData}` containing all the linear constraints extracted from `o`.
"""
function get_all_linear_constraints(o, nvars::Int)
    constraints = Vector{LinearConstraintData}()

    # Process equality constraints.
    eq_ids = MOI.get(o, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}())
    for id in eq_ids
        f = MOI.get(o, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(o, MOI.ConstraintSet(), id)
        rhs, sense = get_rhs_and_type(cons_set)
        b, c0 = build_linear_function_data_sparse(f, nvars)
        push!(constraints, LinearConstraintData(f, b, c0, rhs, sense))
    end

    # Process less-than constraints.
    lt_ids = MOI.get(o, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}())
    for id in lt_ids
        f = MOI.get(o, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(o, MOI.ConstraintSet(), id)
        rhs, sense = get_rhs_and_type(cons_set)
        b, c0 = build_linear_function_data_sparse(f, nvars)
        push!(constraints, LinearConstraintData(f, b, c0, rhs, sense))
    end

    # Process greater-than constraints.
    gt_ids = MOI.get(o, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}}())
    for id in gt_ids
        f = MOI.get(o, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(o, MOI.ConstraintSet(), id)
        rhs, sense = get_rhs_and_type(cons_set)
        b, c0 = build_linear_function_data_sparse(f, nvars)
        push!(constraints, LinearConstraintData(f, b, c0, rhs, sense))
    end

    return constraints
end

"""
    compute_contribution_from_coeff(a, var, lb, ub)

Compute the minimal and maximal contribution for a term with coefficient `a`
associated with variable `var`, given current lower bounds `lb` and upper bounds `ub`.

For a nonnegative coefficient:
  - Minimal contribution: `a * lb[var]`
  - Maximal contribution: `a * ub[var]`

For a negative coefficient:
  - Minimal contribution: `a * ub[var]`
  - Maximal contribution: `a * lb[var]`

Handles infinite bounds by clamping to ±INF_BOUND to avoid NaN/Inf arithmetic issues.

# Arguments
- `a::Float64`: The coefficient.
- `var::Int`: The variable index.
- `lb::Dict{Int,Float64}`: Lower bounds for the variables.
- `ub::Dict{Int,Float64}`: Upper bounds for the variables.

# Returns
A tuple `(cmin, cmax, lb_inf, ub_inf)` with the minimal and maximal contributions,
and booleans indicating if lb or ub were infinite.
"""
function compute_contribution_from_coeff(a::Float64, var::Int, lb::Dict{Int,Float64}, ub::Dict{Int,Float64})
    lb_val = lb[var]
    ub_val = ub[var]

    # Track if bounds are effectively infinite
    lb_inf = lb_val <= -INF_BOUND
    ub_inf = ub_val >= INF_BOUND

    if abs(a) < PROP_TOL_STRICT
        return 0.0, 0.0, false, false
    end

    if a > PROP_TOL_STRICT
        # Positive coefficient: min uses lb, max uses ub
        cmin = lb_inf ? -Inf : a * lb_val
        cmax = ub_inf ? Inf : a * ub_val
    else  # a < -PROP_TOL_STRICT
        # Negative coefficient: min uses ub, max uses lb
        cmin = ub_inf ? -Inf : a * ub_val
        cmax = lb_inf ? Inf : a * lb_val
    end
    return cmin, cmax, lb_inf || (a < 0 && ub_inf), ub_inf || (a < 0 && lb_inf)
end

"""
    create_prop_constraint_data(lin, lb, ub)

Construct a `PropConstraintData` object from a given `LinearConstraintData` (`lin`)
and the current variable bounds `lb` and `ub`. This function uses the sparse
coefficient vector `lin.b` to compute the contributions of each variable in the constraint.

# Arguments
- `lin::LinearConstraintData`: The linear constraint data.
- `lb::Dict{Int,Float64}`: A dictionary of lower bounds for the variables.
- `ub::Dict{Int,Float64}`: A dictionary of upper bounds for the variables.

# Returns
A `PropConstraintData` object with the following fields:
  - `min_activity`: The sum of the minimal contributions (finite part).
  - `max_activity`: The sum of the maximal contributions (finite part).
  - `min_act_inf_cnt`: Count of variables contributing -Inf to min activity.
  - `max_act_inf_cnt`: Count of variables contributing +Inf to max activity.
  - `contributions`: A dictionary mapping variable indices to a tuple `(cmin, cmax)`.
"""
function create_prop_constraint_data(lin::LinearConstraintData,
    lb::Dict{Int,Float64},
    ub::Dict{Int,Float64})
    min_activity = 0.0
    max_activity = 0.0
    min_act_inf_cnt = 0
    max_act_inf_cnt = 0
    contributions = Dict{Int,Tuple{Float64,Float64}}()

    # Iterate over nonzero entries in the sparse coefficient vector.
    # `lin.b.nzind` contains variable indices and `lin.b.nzval` the corresponding coefficients.
    for (var, a) in zip(lin.b.nzind, lin.b.nzval)
        cmin, cmax, min_inf, max_inf = compute_contribution_from_coeff(a, var, lb, ub)
        contributions[var] = (cmin, cmax)
        if min_inf
            min_act_inf_cnt += 1
        else
            min_activity += cmin
        end
        if max_inf
            max_act_inf_cnt += 1
        else
            max_activity += cmax
        end
    end

    return PropConstraintData(lin, min_activity, max_activity, min_act_inf_cnt, max_act_inf_cnt, contributions)
end

"""
    update_constraint_activity!(prop_con, var, lb, ub)

Update the activity of a propagation constraint `prop_con` for a given variable `var`
based on the current bounds in `lb` and `ub`. This function recomputes the minimal and
maximal contribution of `var` and updates the corresponding entries in `prop_con.contributions`
as well as the overall `prop_con.min_activity` and `prop_con.max_activity`.

# Arguments
- `prop_con::PropConstraintData`: The propagation constraint to update.
- `var::Int`: The index of the variable whose bounds have changed.
- `lb::Dict{Int,Float64}`: Dictionary of current lower bounds.
- `ub::Dict{Int,Float64}`: Dictionary of current upper bounds.

# Returns
`true` if the constraint's activity was updated; otherwise, `false`.
"""
function update_constraint_activity!(prop_con::PropConstraintData,
    var::Int,
    lb::Dict{Int,Float64},
    ub::Dict{Int,Float64})
    # Retrieve the coefficient for the variable from the sparse vector.
    a = prop_con.lin.b[var]
    # If the coefficient is effectively zero, nothing to update.
    if abs(a) < PROP_TOL_STRICT
        return false
    end

    old_cmin, old_cmax = prop_con.contributions[var]
    old_min_inf = isinf(old_cmin) && old_cmin < 0
    old_max_inf = isinf(old_cmax) && old_cmax > 0

    new_cmin, new_cmax, new_min_inf, new_max_inf = compute_contribution_from_coeff(a, var, lb, ub)

    updated = false

    # Update min activity
    if old_min_inf != new_min_inf || (!old_min_inf && !new_min_inf && abs(new_cmin - old_cmin) > PROP_TOL_STRICT)
        if old_min_inf && !new_min_inf
            prop_con.min_act_inf_cnt -= 1
            prop_con.min_activity += new_cmin
        elseif !old_min_inf && new_min_inf
            prop_con.min_act_inf_cnt += 1
            prop_con.min_activity -= old_cmin
        else
            prop_con.min_activity += new_cmin - old_cmin
        end
        updated = true
    end

    # Update max activity
    if old_max_inf != new_max_inf || (!old_max_inf && !new_max_inf && abs(new_cmax - old_cmax) > PROP_TOL_STRICT)
        if old_max_inf && !new_max_inf
            prop_con.max_act_inf_cnt -= 1
            prop_con.max_activity += new_cmax
        elseif !old_max_inf && new_max_inf
            prop_con.max_act_inf_cnt += 1
            prop_con.max_activity -= old_cmax
        else
            prop_con.max_activity += new_cmax - old_cmax
        end
        updated = true
    end

    if updated
        prop_con.contributions[var] = (new_cmin, new_cmax)
    end
    return updated
end

"""
    build_var_to_constraints(prop_constraints)

Constructs a mapping from variable indices to a vector of indices of propagation constraints
in which each variable appears. This function iterates over each propagation constraint in
`prop_constraints` and, for each nonzero coefficient (as indicated by the sparse vector),
associates the variable index with the constraint's index.

# Arguments
- `prop_constraints::Vector{PropConstraintData}`: A vector of propagation constraints.

# Returns
A `Dict{Int, Vector{Int}}` mapping each variable index to a vector of constraint indices.
"""
function build_var_to_constraints(prop_constraints::Vector{PropConstraintData})
    var_to_constraints = Dict{Int,Vector{Int}}()
    for (i, prop_con) in enumerate(prop_constraints)
        # Iterate over indices of nonzero entries in the sparse vector.
        for var in prop_con.lin.b.nzind
            push!(get!(var_to_constraints, var, Int[]), i)
        end
    end
    return var_to_constraints
end

"""
    fix_variables_sequence!(prop_constraints, var_to_constraints, lower_bounds, upper_bounds, fix_list, x_ref)

Iteratively fixes the variables in `fix_list` to their corresponding reference values in `x_ref`
by setting both their lower and upper bounds. For each fixed variable, the function updates all
propagation constraints in which the variable appears and then performs incremental propagation
to tighten bounds across the affected constraints.

# Arguments
- `prop_constraints::Vector{PropConstraintData}`: The vector of propagation constraints.
- `var_to_constraints::Dict{Int, Vector{Int}}`: A mapping from variable indices to the indices of
  propagation constraints that mention them.
- `lower_bounds::Dict{Int, Float64}`: A dictionary of current lower bounds for the variables.
- `upper_bounds::Dict{Int, Float64}`: A dictionary of current upper bounds for the variables.
- `fix_list::Vector{Int}`: A vector of variable indices that are to be fixed.
- `x_ref::Vector{Float64}`: A reference vector where the target value for variable `i` is `x_ref[i]`.

# Side Effects
This function updates the `lower_bounds` and `upper_bounds` dictionaries in-place and prints the total
number of constraint updates made during propagation.
# Returns
- The total number of updates performed during propagation.
- A boolean indicating whether an infeasibility was detected during fix and propagate.
"""
function fix_variables_sequence!(prop_constraints::Vector{PropConstraintData},
    var_to_constraints::Dict{Int,Vector{Int}},
    lower_bounds::Dict{Int,Float64},
    upper_bounds::Dict{Int,Float64},
    integer_variables::Vector{Int},
    fix_list::Vector{Int},
    x_ref::Vector{Float64})
    total_updates = 0
    while !isempty(fix_list)
        var = pop!(fix_list)
        # Fix the variable by setting both bounds to the reference value and round if it is an integer variable.
        if var in integer_variables
            # either fix to round(x_ref[var]) or to the closest bound if x_ref[var] is outside the bounds
            val = round(x_ref[var])
            if val < lower_bounds[var]
                lower_bounds[var] = ceil_eps(lower_bounds[var])
                upper_bounds[var] = lower_bounds[var]
            elseif val > upper_bounds[var]
                lower_bounds[var] = floor_eps(upper_bounds[var])
                upper_bounds[var] = lower_bounds[var]
            else
                lower_bounds[var] = val
                upper_bounds[var] = val
            end
        else
            val = x_ref[var]
            if val < lower_bounds[var]
                upper_bounds[var] = lower_bounds[var]
            elseif val > upper_bounds[var]
                lower_bounds[var] = upper_bounds[var]
            else
                lower_bounds[var] = val
                upper_bounds[var] = val
            end
        end
        total_updates += 1
        # Update all constraints that involve this variable.
        relevant_constraints = constraints_for_var(var_to_constraints, var)
        for cid in relevant_constraints
            update_constraint_activity!(prop_constraints[cid], var, lower_bounds, upper_bounds)
        end

        # Propagate the updated bounds across the constraints in which the variable appears.
        number_updates, infeasible = cons_propagate!(prop_constraints, var_to_constraints,
            relevant_constraints, lower_bounds, upper_bounds, integer_variables)
        total_updates += number_updates
        if infeasible
            return 0, true
        end
        length_fix_list = length(fix_list)
        # Remove from fix_list any variable that is already fixed.
        filter!(i -> abs(lower_bounds[i] - upper_bounds[i]) > PROP_TOL_STRICT, fix_list)
        total_updates += length_fix_list - length(fix_list)
    end
    return total_updates, false
end

"""
    propagate_constraint!(prop_con, lb, ub, integer_variables)

Propagates bound changes for a single propagation constraint `prop_con` (of type
`PropConstraintData`), updating the bounds in `lb` and `ub` according to the
constraint's sense (:EqualTo, :LessThan, or :GreaterThan).

For an equality constraint, if only one free (unfixed) variable remains, it is fixed
to the unique value that satisfies the constraint. If more than one free variable remains,
the bounds are narrowed if possible.

# Arguments
- `prop_con::PropConstraintData`: The propagation constraint to update.
- `lb::Dict{Int,Float64}`: Current lower bounds.
- `ub::Dict{Int,Float64}`: Current upper bounds.
- `integer_variables::Vector{Int}`: Indices of integer variables.

# Returns
A `Set{Int}` containing the indices of variables whose bounds were updated, and a boolean
indicating whether the constraint is infeasible.
"""
function propagate_constraint!(prop_con::PropConstraintData,
                               lb::Dict{Int,Float64},
                               ub::Dict{Int,Float64},
                               integer_variables::Vector{Int})
    updated_vars = Set{Int}()
    lin = prop_con.lin
    s = lin.sense    # :EqualTo, :LessThan, or :GreaterThan
    rhs_eff = lin.rhs - lin.f.constant  # effective right-hand side

    if DEBUG_PROPAGATION
        println("Propagating constraint: ", lin.f, " ", s, " ", lin.rhs)
        for var in lin.b.nzind
            println("  Var ", var, ": lb=", lb[var], ", ub=", ub[var])
        end
    end

    # Helper: update bound if the candidate is tighter; update activity contributions accordingly.
    function update_bound!(var::Int, new_lb, new_ub)
        changed = false
        if new_lb > lb[var] + PROP_TOL && new_lb <= ub[var] + PROP_TOL
            lb[var] = new_lb
            changed = true
            if DEBUG_PROPAGATION
                println("  Updated lb for var $var to $new_lb")
            end
        end
        if new_ub < ub[var] - PROP_TOL && new_ub >= lb[var] - PROP_TOL
            ub[var] = new_ub
            changed = true
            if DEBUG_PROPAGATION
                println("  Updated ub for var $var to $new_ub")
            end
        end
        if changed
            push!(updated_vars, var)
            # Update the variable's contribution and overall activity:
            update_constraint_activity!(prop_con, var, lb, ub)
        end
    end

    if s == :EqualTo
        # Feasibility check: need finite bounds on both sides
        if prop_con.min_act_inf_cnt == 0 && prop_con.min_activity > rhs_eff + PROP_TOL
            if DEBUG_PROPAGATION
                println("Infeasible (EqualTo): min_activity=$(prop_con.min_activity) > rhs_eff=$rhs_eff")
            end
            return updated_vars, true
        end
        if prop_con.max_act_inf_cnt == 0 && prop_con.max_activity < rhs_eff - PROP_TOL
            if DEBUG_PROPAGATION
                println("Infeasible (EqualTo): max_activity=$(prop_con.max_activity) < rhs_eff=$rhs_eff")
            end
            return updated_vars, true
        end

        # Identify free (unfixed) variables.
        free_vars = [var for var in lin.b.nzind if abs(lb[var] - ub[var]) > PROP_TOL]

        if length(free_vars) == 1 && prop_con.min_act_inf_cnt == 0
            # Single free variable with finite activity: fix it.
            var = free_vars[1]
            a = lin.b[var]
            current_contrib = prop_con.contributions[var][1]
            others = prop_con.min_activity - current_contrib
            candidate = (rhs_eff - others) / a
            if candidate < lb[var] - PROP_TOL || candidate > ub[var] + PROP_TOL
                return updated_vars, true  # infeasible candidate
            end
            if var in integer_variables && abs(candidate - round(candidate)) > PROP_TOL
                return updated_vars, true  # infeasible candidate
            end
            update_bound!(var, candidate, candidate)
        end
        # TODO: Narrow bounds for multiple free variables

    elseif s == :LessThan
        # For constraints of the form ∑ a_i x_i ≤ rhs_eff,
        # the current minimal activity must not exceed rhs_eff (if finite).
        if prop_con.min_act_inf_cnt == 0 && prop_con.min_activity > rhs_eff + PROP_TOL
            if DEBUG_PROPAGATION
                println("Infeasible (LessThan): min_activity=$(prop_con.min_activity) > rhs_eff=$rhs_eff")
            end
            return updated_vars, true
        end

        # Can only propagate if min activity is finite (no variables with -inf contribution)
        if prop_con.min_act_inf_cnt == 0
            for (var, a) in zip(lin.b.nzind, lin.b.nzval)
                # Skip fixed variables.
                if abs(lb[var] - ub[var]) < PROP_TOL
                    continue
                end
                if abs(a) > PROP_TOL_STRICT
                    contrib_min = prop_con.contributions[var][1]
                    others = prop_con.min_activity - contrib_min
                    beta = rhs_eff - others
                    if a > PROP_TOL_STRICT && a > beta
                        # For a positive coefficient, lowering ub tightens the sum.
                        candidate = lb[var] + beta / a
                        new_ub = var in integer_variables ? floor_eps(candidate) : candidate
                        update_bound!(var, lb[var], new_ub)
                    elseif a < -PROP_TOL_STRICT && -a > beta
                        # For a negative coefficient, raising lb tightens the sum.
                        candidate = ub[var] + beta / a
                        new_lb = var in integer_variables ? ceil_eps(candidate) : candidate
                        update_bound!(var, new_lb, ub[var])
                    end
                end
            end
        end

    elseif s == :GreaterThan
        # For constraints of the form ∑ a_i x_i ≥ rhs_eff,
        # the current maximal activity must not be below rhs_eff (if finite).
        if prop_con.max_act_inf_cnt == 0 && prop_con.max_activity < rhs_eff - PROP_TOL
            if DEBUG_PROPAGATION
                println("Infeasible (GreaterThan): max_activity=$(prop_con.max_activity) < rhs_eff=$rhs_eff")
            end
            return updated_vars, true
        end

        # Can only propagate if max activity is finite (no variables with +inf contribution)
        # FIX: Use max_activity and max contributions, not min_activity
        if prop_con.max_act_inf_cnt == 0
            for (var, a) in zip(lin.b.nzind, lin.b.nzval)
                if abs(lb[var] - ub[var]) < PROP_TOL
                    continue
                end
                if abs(a) > PROP_TOL_STRICT
                    contrib_max = prop_con.contributions[var][2]  # Use max contribution
                    others = prop_con.max_activity - contrib_max   # Use max_activity
                    beta = others - rhs_eff  # surplus = max_activity - rhs

                    # For a positive coefficient: if alpha > beta, raise lb
                    # alpha = a * (ub - lb) for positive a
                    if a > PROP_TOL_STRICT
                        alpha = a * (ub[var] - lb[var])
                        if alpha > beta + PROP_TOL
                            candidate = ub[var] - beta / a
                            new_lb = var in integer_variables ? ceil_eps(candidate) : candidate
                            update_bound!(var, new_lb, ub[var])
                        end
                    elseif a < -PROP_TOL_STRICT
                        # For negative coefficient: if alpha > beta, lower ub
                        # alpha = -a * (ub - lb) for negative a
                        alpha = -a * (ub[var] - lb[var])
                        if alpha > beta + PROP_TOL
                            candidate = lb[var] - beta / a
                            new_ub = var in integer_variables ? floor_eps(candidate) : candidate
                            update_bound!(var, lb[var], new_ub)
                        end
                    end
                end
            end
        end

    else
        error("Unknown constraint sense: $s")
    end

    if DEBUG_PROPAGATION && !isempty(updated_vars)
        println("After propagation for constraint: ", lin.f, " ", s, " ", lin.rhs)
        for var in lin.b.nzind
            println("  Var ", var, ": lb=", lb[var], ", ub=", ub[var])
        end
    end

    return updated_vars, false
end

"""
    cons_propagate!(prop_constraints, var_to_constraints, relevant_constraints, lb, ub)

propagates bound changes across the specified propagation constraints until
no further improvements can be made.

# Arguments
- `prop_constraints::Vector{PropConstraintData}`: A vector of all propagation constraints.
- `var_to_constraints::Dict{Int, Vector{Int}}`: A mapping from variable indices to a vector of
  constraint indices in which they appear.
- `relevant_constraints::Vector{Int}`: The initial set of constraint indices to propagate.
- `lb::Dict{Int, Float64}`: A dictionary of current lower bounds.
- `ub::Dict{Int, Float64}`: A dictionary of current upper bounds.

# Returns
The total number of updates performed during propagation, and a boolean indicating whether
an infeasibility was detected.
"""
function cons_propagate!(prop_constraints::Vector{PropConstraintData},
    var_to_constraints::Dict{Int,Vector{Int}},
    relevant_constraints::Vector{Int},
    lb::Dict{Int,Float64},
    ub::Dict{Int,Float64},
    integer_variables::Vector{Int})
    num_updates = 0
    # Initialize the worklist with the relevant constraint indices.
    worklist = copy(relevant_constraints)
    in_worklist = falses(length(prop_constraints))
    in_worklist[worklist] .= true
    while !isempty(worklist)
        i = pop!(worklist)
        in_worklist[i] = false
        local_pc = prop_constraints[i]
        # Propagate the current constraint.
        updated_vars, infeasible = propagate_constraint!(local_pc, lb, ub, integer_variables)
        if infeasible
            return 0, true
        end
        num_updates += length(updated_vars)
        # For each updated variable, update all constraints where it appears.
        for var in updated_vars
            for j in constraints_for_var(var_to_constraints, var)
                if update_constraint_activity!(prop_constraints[j], var, lb, ub)
                    if !in_worklist[j]
                        push!(worklist, j)
                        in_worklist[j] = true
                    end
                end
            end
        end
    end
    return num_updates, false
end

"""
    apply_prop(problem_data, fix_list, x_ref)

Extracts linear constraints from `problem_data`, builds the corresponding propagation
constraints and a mapping from variables to constraints, fixes the variables specified in
`fix_list` (using a reference vector of zeros), and applies incremental propagation.

# Arguments
- `problem_data`: A structure containing the original model, optimizer, and global variable bounds.
- `fix_list::Vector{Int}`: A vector of variable indices to fix.
- `x_ref::Vector{Float64}`: A reference vector used to fix the variables.

# Returns
A tuple `(lower_bounds, upper_bounds)` with the updated variable bounds.
"""
function apply_prop(problem_data, fix_list, x_ref)
    # Determine the number of variables.
    num_vars = num_variables(problem_data.optimizer)
    # Extract the indices of integer variables.
    integer_variables = get_integer_variables_indices(problem_data.optimizer)
    # Extract linear constraints from the optimizer.
    linear_constraints = get_all_linear_constraints(problem_data.optimizer, num_vars)

    # Copy global variable bounds.
    lower_bounds = copy(problem_data.global_variable_bounds.lower_bounds)
    upper_bounds = copy(problem_data.global_variable_bounds.upper_bounds)

    # Create propagation constraints using the sparse representation.
    prop_constraints = [create_prop_constraint_data(lin, lower_bounds, upper_bounds)
                        for lin in linear_constraints]

    # Build the mapping from variable indices to constraint indices.
    var_to_constraints = build_var_to_constraints(prop_constraints)

    # copy the list of variables to fix
    to_fix = copy(fix_list)

    # Fix the variables in `fix_list` and propagate the bound changes.
    total_updates, infeasible = fix_variables_sequence!(prop_constraints, var_to_constraints, lower_bounds, upper_bounds, integer_variables, to_fix, x_ref)
    if DEBUG_PROPAGATION
        println("Total number of updates: ", total_updates)
        println("Infeasible: ", infeasible)
    end
    @assert infeasible || (total_updates >= length(fix_list)) "Didn't fix enough variables: $total_updates < $(length(fix_list))"

    return lower_bounds, upper_bounds, infeasible, total_updates
end

"""
    propagate_all!(prop_constraints, var_to_constraints, lower_bounds, upper_bounds)

Propagates bound changes over all constraints until no further improvements (bound reductions) can be made.
It initializes the propagation worklist with all constraint indices and calls the `cons_propagate!` function.
Returns the total number of updates performed.
"""
function propagate_all(optimizer, lbs, ubs, num_vars)

    # Extract linear constraints from the optimizer.
    linear_constraints = get_all_linear_constraints(optimizer, num_vars)

    # Extract the indices of integer variables.
    integer_variables = get_integer_variables_indices(optimizer)

    # Copy global variable bounds.
    lower_bounds = copy(lbs)
    upper_bounds = copy(ubs)

    # Create propagation constraints using the sparse representation.
    prop_constraints = [create_prop_constraint_data(lin, lower_bounds, upper_bounds) for lin in linear_constraints]

    # Build the mapping from variable indices to constraint indices.
    var_to_constraints = build_var_to_constraints(prop_constraints)

    # Create a worklist with all propagation constraints.
    all_constraints = collect(1:length(prop_constraints))
    total_updates, infeasible = cons_propagate!(prop_constraints, var_to_constraints, all_constraints, lower_bounds, upper_bounds, integer_variables)
    if DEBUG_PROPAGATION
        println("Propagation result feasible: ", !infeasible)
        println("Total number of updates: ", total_updates)
    end
    for v in keys(lower_bounds)
        if lower_bounds[v] > upper_bounds[v]
            @assert false "Infeasible bounds!"
        end
    end
    return lower_bounds, upper_bounds, infeasible, total_updates
end

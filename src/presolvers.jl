function transform_problem_mccormick!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model, bigM)
    optimizer = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    MOI.copy_to(optimizer, linear_model)

    norig_vars = num_variables(linear_model)

    v_bounds = build_global_bounds(optimizer, bigM)
    lbs = v_bounds.lower_bounds
    ubs = v_bounds.upper_bounds

    # here we have ones if xi * xj is in the model
    # and zeros if xi * xj does not appear in the model
    interaction_matrix = create_sparse_quadratic_interaction_matrix(num_variables(linear_model), (spzeros(0, 0), Float64[], 0.0), quadratic_lowerthan, quadratic_equalto)

    println("Number of variables in linear model: ", MOI.get(linear_model, MOI.NumberOfVariables()))

    # create new variables continuous variable in a dictionary (i,j) -> z_{i,j}
    z_vars = Dict{Tuple{Int,Int},MOI.VariableIndex}()
    # consider only upper diagonal
    for (i, j) in zip(findnz(interaction_matrix)...)
        if i > j
            continue
        end
        ordered_i, ordered_j = sort([i, j])
        lb_z_var = min(lbs[i] * lbs[j], lbs[i] * ubs[j], ubs[i] * lbs[j], ubs[i] * ubs[j])
        ub_z_var = max(lbs[i] * lbs[j], lbs[i] * ubs[j], ubs[i] * lbs[j], ubs[i] * ubs[j])
        z_bounds = (MOI.GreaterThan(lb_z_var), MOI.LessThan(ub_z_var))
        z_var, (_, _) = MOI.add_constrained_variable(original_model, z_bounds)
        z_var, (_, _) = MOI.add_constrained_variable(linear_model, z_bounds)
        z_vars[(ordered_i, ordered_j)] = z_var
    end
    number_added_vars = length(z_vars)
    println("Number of added variables in the mccormick model: ", number_added_vars)
    @assert MOI.get(linear_model, MOI.NumberOfVariables()) == norig_vars + number_added_vars

    # For every quadratic constraint add a linear constraint where we replace all xi * xj by zij
    # E.g. xi * xj + xk * xl = 3 -> zij + zkl = 3.
    for (c_idx, (Q, b, c0)) in quadratic_equalto
        new_b = copy(b)
        # append zeros for the new variables
        append!(new_b, zeros(length(z_vars)))
        new_c0 = c0
        for (i, j) in zip(findnz(Q)...)
            ordered_i, ordered_j = sort([i, j])
            z_var = z_vars[(ordered_i, ordered_j)]
            # get index of z_var
            index_z_var = z_var.value
            new_b[index_z_var] = 2.0 * Q[i, j]
        end
        # Create a vector of scalar affine terms for each variable in x.
        affine_terms = [MOI.ScalarAffineTerm(new_b[i], MOI.VariableIndex(i)) for i in eachindex(new_b)]

        # Construct the scalar affine function representing new_b * x
        aff_func = MOI.ScalarAffineFunction(affine_terms, 0.0)

        # Create the constraint new_b*x == -new_c0.
        MOI.add_constraint(linear_model, aff_func, MOI.EqualTo(-new_c0))
    end
    for (c_idx, (Q, b, c0)) in quadratic_lowerthan
        new_b = copy(b)
        # append zeros for the new variables
        append!(new_b, zeros(length(z_vars)))
        new_c0 = c0
        for (i, j) in zip(findnz(Q)...)
            ordered_i, ordered_j = sort([i, j])
            z_var = z_vars[(ordered_i, ordered_j)]
            # get index of z_var
            index_z_var = z_var.value
            new_b[index_z_var] = 2.0 * Q[i, j]
        end
        # Create a vector of scalar affine terms for each variable in x.
        affine_terms = [MOI.ScalarAffineTerm(new_b[i], MOI.VariableIndex(i)) for i in eachindex(new_b)]

        # Construct the scalar affine function representing new_b * x.
        aff_func = MOI.ScalarAffineFunction(affine_terms, 0.0)

        # Create the constraint new_b*x <= -new_c0 <= 0.
        MOI.add_constraint(linear_model, aff_func, MOI.LessThan(-new_c0))
    end

    num_added_constraints = 0
    # For every new variable add the four mccormick inequalities
    # zij >= lb(xi) * xj + xi * lb(xj) - lb(xi) * lb(xj)
    # zij >= ub(xi) * xj + xi * ub(xj) - ub(xi) * ub(xj)
    # zij <= ub(xi) * xj + xi * lb(xj) - ub(xi) * lb(xj)
    # zij <= xi * ub(xj) + lb(xi) * xj - lb(xi) * ub(xj)
    # and only if the lb and ub are finite
    for ((i, j), z_var) in z_vars
        lb_i = lbs[i]
        ub_i = ubs[i]
        lb_j = lbs[j]
        ub_j = ubs[j]
        # if any of the bound products below is larger than bigM we do not add the constraint
        if abs(lb_i * lb_j) >= bigM || abs(lb_i * ub_j) >= bigM || abs(ub_i * lb_j) >= bigM || abs(ub_i * ub_j) >= bigM
            continue
        end
        # zij - lb(xi) * xj - xi * lb(xj) >= - lb(xi) * lb(xj)
        terms = [
            MOI.ScalarAffineTerm(1.0, z_var),
            MOI.ScalarAffineTerm(-lb_i, MOI.VariableIndex(j)),
            MOI.ScalarAffineTerm(-lb_j, MOI.VariableIndex(i))
        ]
        MOI.add_constraint(linear_model, MOI.ScalarAffineFunction(terms, 0.0), MOI.GreaterThan(-lb_i * lb_j))

        # zij - ub(xi) * xj - xi * ub(xj) >= - ub(xi) * ub(xj)
        terms = [
            MOI.ScalarAffineTerm(1.0, z_var),
            MOI.ScalarAffineTerm(-ub_i, MOI.VariableIndex(j)),
            MOI.ScalarAffineTerm(-ub_j, MOI.VariableIndex(i))
        ]
        MOI.add_constraint(linear_model, MOI.ScalarAffineFunction(terms, 0.0), MOI.GreaterThan(-ub_i * ub_j))

        # zij - ub(xi) * xj - xi * lb(xj) <= - ub(xi) * lb(xj)
        terms = [
            MOI.ScalarAffineTerm(1.0, z_var),
            MOI.ScalarAffineTerm(-ub_i, MOI.VariableIndex(j)),
            MOI.ScalarAffineTerm(-lb_j, MOI.VariableIndex(i))
        ]
        MOI.add_constraint(linear_model, MOI.ScalarAffineFunction(terms, 0.0), MOI.LessThan(-ub_i * lb_j))

        # zij - lb(xi) * xj - xi * ub(xj) <= - lb(xi) * ub(xj)
        terms = [
            MOI.ScalarAffineTerm(1.0, z_var),
            MOI.ScalarAffineTerm(-lb_i, MOI.VariableIndex(j)),
            MOI.ScalarAffineTerm(-ub_j, MOI.VariableIndex(i))
        ]
        MOI.add_constraint(linear_model, MOI.ScalarAffineFunction(terms, 0.0), MOI.LessThan(-lb_i * ub_j))
        num_added_constraints += 4
    end
    println("Number of added constraints in the mccormick model: ", num_added_constraints)
    Q_obj, b_obj, c0_obj = quad_objective_data
    b_obj_mccormick = vcat(b_obj, zeros(number_added_vars))
    if length(Q_obj) > 0
        Q_obj_mccormick = vcat(hcat(Q_obj, zeros(size(Q_obj, 1), number_added_vars)), zeros(norig_vars + number_added_vars, number_added_vars)')
    else
        Q_obj_mccormick = copy(Q_obj)
    end
    # update the objective data
    quad_objective_data_mccormick = Q_obj_mccormick, b_obj_mccormick, c0_obj

    @assert MOI.get(linear_model, MOI.NumberOfVariables()) == norig_vars + number_added_vars

    # More efficient code for resizing constraint matrices
    quadratic_lowerthan_mccormick = typeof(quadratic_lowerthan)()
    quadratic_equalto_mccormick = typeof(quadratic_equalto)()
    total_vars = norig_vars + number_added_vars
    for (data_dict, new_dict) in zip([quadratic_lowerthan, quadratic_equalto], [quadratic_lowerthan_mccormick, quadratic_equalto_mccormick])
        for (c_idx, (Q, b, c0)) in data_dict
            # Create new extended b vector
            new_b = zeros(total_vars)
            @views new_b[1:norig_vars] .= b

            # Create new extended sparse Q matrix more efficiently
            if isa(Q, SparseMatrixCSC)
                rows, cols, vals = findnz(Q)
                new_Q = sparse(rows, cols, vals, total_vars, total_vars)
            else
                # Handle the case where Q might be dense
                new_Q = spzeros(total_vars, total_vars)
                @views new_Q[1:norig_vars, 1:norig_vars] .= Q
            end

            @assert size(new_Q) == (total_vars, total_vars)
            @assert length(new_b) == total_vars

            new_dict[c_idx] = (new_Q, new_b, c0)
        end
    end
    return quadratic_lowerthan_mccormick, quadratic_equalto_mccormick, quad_objective_data_mccormick, number_added_vars, num_added_constraints
end

function mccormick_relaxation!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model, bigM)
    # changes all problem data in place
    quadratic_lowerthan_mccormick, quadratic_equalto_mccormick, quad_objective_data_mccormick, number_added_vars, number_added_constraints =
                                                     transform_problem_mccormick!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model, bigM)
    function postsolve_mccormick(solution::Vector{Float64})
        new_solution = solution[1:(end-number_added_vars)]
        return new_solution
    end
    return number_added_vars, number_added_constraints, quadratic_lowerthan_mccormick, quadratic_equalto_mccormick, quad_objective_data_mccormick, postsolve_mccormick
end

function presolve_bilinear!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model)
    # changes all problem data in place
    quad_objective_data, number_added_vars = transform_problem_complementarities!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model)
    function postsolve_bilinear(solution::Vector{Float64})
        new_solution = solution[1:(end-number_added_vars)]
        return new_solution
    end
    return quad_objective_data, postsolve_bilinear
end

function transform_problem_complementarities!(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model)
    Q_obj, b_obj, c0_obj = quad_objective_data
    bilinear_complementarities = []
    for (c_idx, (Q, b, c0)) in quadratic_equalto
        if SparseArrays.nnz(Q) == 2 && norm(b) ≈ 0 && abs(c0) ≈ 0
            # we extract the two indices from the row
            row_indices, _, _ = SparseArrays.findnz(Q)
            @assert length(row_indices) == 2
            x1, x2 = row_indices
            push!(bilinear_complementarities, (c_idx, x1, x2))
        end
    end
    @info "Reformulating $(length(bilinear_complementarities)) complementarities"
    norig_vars = length(b_obj)
    fixed_vars = Int[]
    # add auxiliary binary variables and remove the quadratic from both models
    for (c_idx, x1, x2) in bilinear_complementarities
        x1_upperbound = if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x1))
            MOI.get(linear_model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x1)).upper
        else
            Inf
        end
        x1_lowerbound = if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x1))
            MOI.get(linear_model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x1)).lower
        else
            -Inf
        end
        constraint_set1 = if x1_lowerbound >= 0
            MOI.LessThan(0.0)
        elseif x1_upperbound <= 0
            MOI.GreaterThan(0.0)
        else
            MOI.EqualTo(0.0)
        end
        x2_upperbound = if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2))
            MOI.get(linear_model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2)).upper
        else
            Inf
        end
        x2_lowerbound = if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2))
            MOI.get(linear_model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2)).lower
        else
            -Inf
        end
        constraint_set2 = if x2_lowerbound >= 0
            MOI.LessThan(0.0)
        elseif x2_upperbound <= 0
            MOI.GreaterThan(0.0)
        else
            MOI.EqualTo(0.0)
        end
        # verify whether we are in a trivial case
        # one variable bounded away from zero -> other set to zero
        if x1_lowerbound > 1e-6 || x1_upperbound < -1e-6
            @assert x2_lowerbound <= 0
            @assert x2_upperbound >= 0
            for model in (original_model, linear_model)
                if MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2))
                    MOI.set(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2), MOI.LessThan{Float64}(0.0))
                else
                    MOI.add_constraint(model, MOI.VariableIndex(x2), MOI.LessThan(0.0))
                end
                if MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2))
                    MOI.set(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2), MOI.GreaterThan{Float64}(0.0))
                else
                    MOI.add_constraint(model, MOI.VariableIndex(x2), MOI.GreaterThan(0.0))
                end
            end
            push!(fixed_vars, x2)
        elseif x2_lowerbound > 1e-6 || x2_upperbound < -1e-6
            @assert x1_lowerbound <= 0
            @assert x1_upperbound >= 0
            push!(fixed_vars, x1)
        end
        if x1 in fixed_vars || x2 in fixed_vars
            continue
        end
        for model in (original_model, linear_model)
            z_new, _ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            vec_func_var1 = MOI.VectorAffineFunction(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, z_new)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(x1)))
                ],
                [0.0, 0.0]
            )
            MOI.add_constraint(model, vec_func_var1, MOI.Indicator{MOI.ACTIVATE_ON_ONE}(constraint_set1))
            vec_func_var2 = MOI.VectorAffineFunction(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, z_new)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(x2)))
                ],
                [0.0, 0.0]
            )
            MOI.add_constraint(model, vec_func_var2, MOI.Indicator{MOI.ACTIVATE_ON_ZERO}(constraint_set2))
            # TODO add optional valid inequalities for convex hull of the complementarity constraint
            # set variable name to be able to delete it at the end
            MOI.set(model, MOI.VariableName(), z_new, "auxvar_bilinear_$(z_new.value)")
            # delete equality constraint which has been replaced
        end
        delete!(quadratic_equalto, c_idx)
    end
    nbins_new = length(bilinear_complementarities) - length(fixed_vars)
    @assert MOI.get(original_model, MOI.NumberOfVariables()) == norig_vars + nbins_new
    @assert MOI.get(linear_model, MOI.NumberOfVariables()) == norig_vars + nbins_new
    # we keep variables fixed to zero for simplicity in the objective and constraints
    # first edit the objective data to account for the new binary variables and the removed fixed variables
    append!(b_obj, zeros(nbins_new))
    if length(Q_obj) > 0
        Q_obj = vcat(hcat(Q_obj, zeros(size(Q_obj, 1), nbins_new)), zeros(norig_vars + nbins_new, nbins_new)')
    end
    # perform the same for all constraints
    for data_dict in (quadratic_equalto, quadratic_lowerthan)
        for (c_idx, (Q, b, c0)) in data_dict
            b = vcat(b, zeros(nbins_new))
            Q = vcat(hcat(Q, zeros(size(Q, 1), nbins_new)), zeros(norig_vars + nbins_new, nbins_new)')
            @assert size(Q) == (norig_vars + nbins_new, norig_vars + nbins_new)
            @assert size(b) == (norig_vars + nbins_new,)
            data_dict[c_idx] = (Q, b, c0)
        end
    end
    quad_objective_data = Q_obj, b_obj, c0_obj
    return quad_objective_data, nbins_new
end

function presolve_deperspective(linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model; use_logbarrier = false, use_indicator = true)
    # changes all problem data in place
    Q_obj, c_obj, c0_obj = quad_objective_data
    perspective_cons_idx = []
    epigraph_indices_info = []
    # first search all perspective candidates
    for (key, (Q_con, b_con, c0_con)) in quadratic_lowerthan
        if SparseArrays.nnz(Q_con) != 3 || norm(b_con) > 1e-6 || abs(c0_con) > 1e-6
            continue
        end
        # only one element on the diagonal
        if SparseArrays.nnz(diag(Q_con)) != 1
            continue
        end
        cols_idx, rows_idx, vals = findnz(Q_con)
        # there should be one left-hand side and two right-hand-side entries
        if count(<(0), vals) != 2 || count(>(0), vals) != 1
            continue
        end
        diag_idx = findfirst(>(0), vals)
        # the left-hand side should be a diagonal term xᵢ²
        if cols_idx[diag_idx] != rows_idx[diag_idx]
            continue
        end
        diag_coeff = vals[diag_idx]
        varidx_squared = rows_idx[diag_idx]
        idx_offdiag = diag_idx == 1 ? 2 : 1
        binvar_idx = -1
        epi_idx = -1
        if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(rows_idx[idx_offdiag]))
            binvar_idx = rows_idx[idx_offdiag]
            epi_idx = cols_idx[idx_offdiag]
        elseif MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(cols_idx[idx_offdiag]))
            binvar_idx = cols_idx[idx_offdiag]
            epi_idx = rows_idx[idx_offdiag]
        else # none of the two variables is a binary -> not a perspective constraint -> abort
            continue
        end
        # the epigraph should not be integer or binary
        if MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(epi_idx)) || MOI.is_valid(linear_model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(epi_idx))
            continue
        end
        binvar_coeff = vals[idx_offdiag]
        # the epigraph variable should not appear in a quadratic term in the objective
        if length(Q_obj) > 0 && norm(Q_obj[epi_idx, :]) > 0
            continue
        end
        # the epigraph variable should not appear in other quadratic constraints
        for (key2, (Q_con2, b_con2, _)) in quadratic_lowerthan
            if key != key2
                if abs(b_con2[epi_idx]) > 0 || norm(Q_con2[epi_idx, :]) > 0
                    continue
                end
            end
        end
        # nor should it appear in the linear constraints
        epigraph_in_multiple_constraints = false
        for (F, S) in MOI.get(linear_model, MOI.ListOfConstraintTypesPresent())
            if F <: MOI.ScalarAffineFunction # TODO check VectorAffineFunction too?
                for c_idx in MOI.get(linear_model, MOI.ListOfConstraintIndices{F,S}())
                    func_test = MOI.get(linear_model, MOI.ConstraintFunction(), c_idx)
                    if any(term -> term.variable == MOI.VariableIndex(epi_idx), func_test.terms)
                        @info "Epigraph variable appears in constraint of type $F, $S, ignoring"
                        epigraph_in_multiple_constraints = true
                        break
                    end
                end
                if epigraph_in_multiple_constraints
                    break
                end
            end
        end
        if epigraph_in_multiple_constraints
            continue
        end
        # last check, the epigraph should not be zero in the objective
        # if this happens and we didn't miss anything else, the whole perspective constraint is useless
        # we consider for now it shouldn't happen
        epi_objective_coeff = c_obj[epi_idx]
        if epi_objective_coeff == 0
            @warn "epigraph $epi_idx objective coefficient is zero, something odd in the model, could be presolved?"
            continue
        end
        # we have a constraint coming from a perspective of a quadratic expression in the objective
        push!(perspective_cons_idx, key)
        # in the final solution, epigraph = a * xᵢ² = expression(x) in the objective
        # we store the two indices and the coefficient to be able to postsolve
        # we transform the linear term in the objective b_i * epigraph into b_i * xᵢ²
        push!(
            epigraph_indices_info,
            (epi_idx = epi_idx, binvar_idx = binvar_idx, varidx_squared = varidx_squared, coefficient = -diag_coeff / binvar_coeff * 0.5)
        )
    end
    if !isempty(perspective_cons_idx)
        if length(Q_obj) == 0
            Q_obj = spzeros(length(c_obj), length(c_obj))
        end
        @info "Transforming $(length(perspective_cons_idx)) perspective constraints"
        for idx in eachindex(perspective_cons_idx)
            delete!(quadratic_lowerthan, perspective_cons_idx[idx])
            nt = epigraph_indices_info[idx]
            # transfer the linear epigraph to a diagonal quadratic
            Q_obj[nt.varidx_squared, nt.varidx_squared] = 2.0 * c_obj[nt.epi_idx]
            c_obj[nt.epi_idx] = 0.0
            for model in (linear_model, original_model)
                # we first compute an upper bound for the variable either from itself, from the epigraph, or we make up one
                continuous_upper_bound = 1e5
                if MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(nt.varidx_squared))
                    continuous_upper_bound = MOI.get(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(nt.varidx_squared)).upper
                elseif MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(nt.epi_idx))
                    continuous_upper_bound = MOI.get(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(nt.epi_idx)).upper
                end

                # x - M * b <= 0
                MOI.add_constraint(model, 1.0 * MOI.VariableIndex(nt.varidx_squared) - continuous_upper_bound * MOI.VariableIndex(nt.binvar_idx), MOI.LessThan(0.0))
                if use_indicator
                    # additional indicator constraint: b = 0 -> x <= 0
                    vec_func_var = MOI.VectorAffineFunction(
                        [
                            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(nt.binvar_idx))),
                            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(nt.varidx_squared)))
                        ],
                        [0.0, 0.0]
                    )
                    MOI.add_constraint(model, vec_func_var, MOI.Indicator{MOI.ACTIVATE_ON_ZERO}(MOI.LessThan(0.0)))
                end
                # we fix the epigraph bounds if none exists
                if !MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(nt.epi_idx))
                    MOI.add_constraint(model, MOI.VariableIndex(nt.epi_idx), MOI.LessThan(continuous_upper_bound))
                end
                if !MOI.is_valid(model, MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(nt.epi_idx))
                    MOI.add_constraint(model, MOI.VariableIndex(nt.epi_idx), MOI.GreaterThan(0.0))
                end
            end
        end
    end

    function postsolve_deperspective(solution::Vector{Float64})
        if isempty(perspective_cons_idx)
            return solution
        end
        new_solution = copy(solution)
        # we need to set the epigraph variable to its correct value
        for nt in epigraph_indices_info
            new_solution[nt.epi_idx] = new_solution[nt.varidx_squared]^2 * nt.coefficient
        end
        return new_solution
    end
    # repack the modified objective data
    quad_objective_data = Q_obj, c_obj, c0_obj
    return linear_model, quadratic_lowerthan, quadratic_equalto, quad_objective_data, original_model, postsolve_deperspective
end

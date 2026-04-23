"""
Creates a MOI model from a file.
"""
function create_moi_model(file_path)
    model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    try
        MOI.read_from_file(model, file_path)
    catch
        error("Could not read the file: $file_path")
    end
    return model
end

# CONVENTION: the function is q(x) = 1/2 * x^T Q x + b^T x + c0
function build_quadratic_function_data(f::MOI.ScalarQuadraticFunction, nvars::Int)
    c0 = f.constant
    b = zeros(nvars)
    for term in f.affine_terms
        b[term.variable.value] = term.coefficient
    end
    Q = spzeros(nvars, nvars)
    for term in f.quadratic_terms
        Q[term.variable_1.value, term.variable_2.value] += term.coefficient
        Q[term.variable_2.value, term.variable_1.value] += term.coefficient
        if term.variable_1 == term.variable_2
            Q[term.variable_1.value, term.variable_2.value] /= 2
        end
    end
    return Q, b, c0
end

function build_linear_function_data(f::MOI.ScalarAffineFunction, nvars::Int)
    c0 = f.constant
    b = zeros(nvars)
    for term in f.terms
        b[term.variable.value] = term.coefficient
    end
    return b, c0
end

# return a list of the binary indices, then a list of the integer indices
function get_binary_integer(model)
    list_binary_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}())
    list_integer_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.Integer}())
    return getproperty.(list_binary_cons, :value), getproperty.(list_integer_cons, :value)
end

"""
Copies the original `src_model` into a destination `dst_model`.
The function removes all quadratic constraints and the quadratic objective.
Quadratic constraints are standardized to q(x) <= 0 and q(x) == 0.
Return:
- The destination model `dst_model`
- The vector of quadratic constraints of the form q(x) <= 0
- The vector of equality quadratic constraints of the form q(x) == 0
- The tuple with the objective in minimization form q(x), OR `nothing` if the objective was not quadratic
"""
function copy_mip_new_model(src_model, dst_model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()))
    # safety first: destination should be empty
    if !MOI.is_empty(dst_model)
        @warn("destination model not empty")
        MOI.empty!(dst_model)
    end
    quadratic_lowerthan = Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}},Tuple}()
    quadratic_equalto = Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}},Tuple}()
    custom_filter = build_model_filter_mip(quadratic_lowerthan, quadratic_equalto, src_model)
    filtered_src = MOI.Utilities.ModelFilter(custom_filter, src_model)

    MOI.copy_to(dst_model, filtered_src)
    # if quadratic objective, extract it
    FT = MOI.get(dst_model, MOI.ObjectiveFunctionType())
    f = MOI.get(dst_model, MOI.ObjectiveFunction{FT}())
    nvars = MOI.get(dst_model, MOI.NumberOfVariables())
    if FT <: MOI.ScalarQuadraticFunction
        Q, b, c0 = build_quadratic_function_data(f, nvars)
        # if linear objective, extract it
    elseif FT <: MOI.ScalarAffineFunction
        Q = spzeros(0, 0)
        b, c0 = build_linear_function_data(f, nvars)
    else
        error("Unsupported objective function, neither linear nor quadratic")
    end
    # adapt to min/max sense
    if MOI.get(dst_model, MOI.ObjectiveSense()) == MOI.MAX_SENSE
        Q.nzval .*= -1
        b .*= -1
    end

    # replace quadratic objective with a dummy linear one
    f0 = zero(MOI.ScalarAffineFunction{Float64})
    MOI.set(dst_model, MOI.ObjectiveFunction{typeof(f0)}(), f0)

    return dst_model, quadratic_lowerthan, quadratic_equalto, (Q, b, c0)
end

function build_model_filter_mip(quadratic_lowerthan::Dict, quadratic_equalto::Dict, model)
    custom_filter(::Any) = true
    nvars = MOI.get(model, MOI.NumberOfVariables())
    function custom_filter(ci::MOI.ConstraintIndex{<:MOI.ScalarQuadraticFunction,S}) where {S}
        f = MOI.get(model, MOI.ConstraintFunction(), ci)
        q_data = build_quadratic_function_data(f, nvars)
        cons_set = MOI.get(model, MOI.ConstraintSet(), ci)
        Q, b, c0 = q_data
        is_equality = false
        if cons_set isa MOI.LessThan
            # q(x) = 1/2 * x^T Q x + b^T x + c0 <= c
            # <=>
            # 1/2 * x^T Q x + b^T x + c0 - c <= 0
            c0 -= cons_set.upper
        elseif cons_set isa MOI.GreaterThan
            # q(x) = 1/2 * x^T Q x + b^T x + c0 >= c
            # <=>
            # -1/2 * x^T Q x - b^T x - c0 + c <= 0
            Q.nzval .*= -1
            b .*= -1
            c0 = cons_set.lower - c0
        elseif cons_set isa MOI.EqualTo
            # normalize to:
            # 1/2 * x^T Q x + b^T x + c0 - c == 0
            c0 -= cons_set.value
            is_equality = true
        else
            error("Unknown set: $(cons_set)")
        end
        new_q_data = (Q, b, c0)
        if is_equality
            push!(quadratic_equalto, ci => new_q_data)
        else
            push!(quadratic_lowerthan, ci => new_q_data)
        end
        return false
    end
    return custom_filter
end

# Function to create the interaction matrix as a sparse array
function create_sparse_quadratic_interaction_matrix(n, quad_objective_data, quadratic_lowerthan, quadratic_equalto)
    interaction_matrix = spzeros(Int, n, n)

    Q_obj, c_obj, c0_obj = quad_objective_data

    if !isempty(Q_obj)
        for (i, j) in zip(findnz(Q_obj)...)
            interaction_matrix[i, j] = 1
            interaction_matrix[j, i] = 1
        end
    end

    for (key, val) in quadratic_lowerthan
        Q_con, b_con, c0_con = val
        for (i, j) in zip(findnz(Q_con)...)
            interaction_matrix[i, j] = 1
            interaction_matrix[j, i] = 1
        end
    end

    for (key, val) in quadratic_equalto
        Q_con, b_con, c0_con = val
        for (i, j) in zip(findnz(Q_con)...)
            interaction_matrix[i, j] = 1
            interaction_matrix[j, i] = 1
        end
    end

    return interaction_matrix
end

"""
Solve a MIP to compute the min cover.
"""
function min_vertex_cover(matrix::SparseMatrixCSC{Int,Int}, solver::String)
    n = size(matrix, 1)
    if solver == "Gurobi"
        optimizer = Gurobi.Optimizer()
    else
        error("Unknown solver: $solver")
    end
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.TimeLimitSec(), 5.0)

    c = ones(n)
    # min 1'x
    # s.t. x_i + x_j >= 1 for eachedge (i, j)
    # x binary
    x = MOI.add_variables(optimizer, length(c))
    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
    )
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # For each edge (i, j) add the constraint x_i + x_j >= 1
    for i in 1:n
        row_indices, _ = findnz(matrix[i, :])
        for j in row_indices
            if i <= j
                MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i]), MOI.ScalarAffineTerm(1.0, x[j])], 0.0), MOI.GreaterThan(1.0))
            end
        end
    end

    # add integrality constraints
    for x_i in x
        MOI.add_constraint(optimizer, x_i, MOI.ZeroOne())
    end

    MOI.optimize!(optimizer)

    if MOI.get(optimizer, MOI.ResultCount()) == 0
        return []
    end

    # return only the indices of the selected vertices
    return findall(MOI.get(optimizer, MOI.VariablePrimal(), x) .> 0.5)
end

"""
Generate multiple vertex covers for the undercover heuristic.
Return a list of different covers.
"""
function diverse_vertex_covers(matrix::SparseMatrixCSC{Int,Int}, solver::String;
                              num_covers=1, max_overlap_ratio=0.7, time_limit=5.0)
    n = size(matrix, 1)
    covers = Vector{Vector{Int}}()

    # First get the minimum vertex cover as the base solution
    min_cover = min_vertex_cover(matrix, solver)
    if !isempty(min_cover)
        push!(covers, min_cover)
    else
        return covers
    end

    if solver == "Gurobi"
        optimizer = Gurobi.Optimizer()
    else
        error("Unknown solver: $solver")
    end

    # Settings for the optimizer
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.TimeLimitSec(), time_limit)

    for i in 1:num_covers-1
        MOI.empty!(optimizer)

        # Set up optimization problem
        x = MOI.add_variables(optimizer, n)

        # Add binary constraints
        for x_i in x
            MOI.add_constraint(optimizer, x_i, MOI.ZeroOne())
        end

        # Add edge covering constraints
        for i in 1:n
            row_indices, _ = findnz(matrix[i, :])
            for j in row_indices
                if i <= j  # Only add once per edge
                    MOI.add_constraint(optimizer,
                        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i]), MOI.ScalarAffineTerm(1.0, x[j])], 0.0),
                        MOI.GreaterThan(1.0))
                end
            end
        end

        # Add diversity constraints based on previous covers
        for prev_cover in covers
            # Limit the overlap with previous covers
            max_overlap = ceil(Int, length(prev_cover) * max_overlap_ratio)

            # Create a constraint to limit overlap with this previous cover
            overlap_terms = [MOI.ScalarAffineTerm(1.0, x[idx]) for idx in prev_cover]
            if !isempty(overlap_terms)
                MOI.add_constraint(optimizer,
                    MOI.ScalarAffineFunction(overlap_terms, 0.0),
                    MOI.LessThan(Float64(max_overlap)))
            end
        end

        # Objective: minimize size but with some randomness
        weights = rand(n) .+ 0.5  # Random weights between 0.5 and 1.5
        MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(weights, x), 0.0))
        MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

        # Solve
        MOI.optimize!(optimizer)

        if MOI.get(optimizer, MOI.ResultCount()) > 0
            # Get the new cover
            new_cover = findall(MOI.get(optimizer, MOI.VariablePrimal(), x) .> 0.5)
            if !isempty(new_cover)
                push!(covers, new_cover)
            end
        end
    end

    return covers
end

# check if the coverset is only on integer variables
function is_cover_set_only_integer(cover_set, integer_variables)
    for col in cover_set
        if !(col in integer_variables)
            return false
        end
    end
    return true
end

function extract_perspective_cones(model, quadratic_lower_than, quad_objective_data)
    Q_obj, c_obj, c0_obj = quad_objective_data
    perspective_cons_idx = Int[]
    # first search all perspective candidates
    for (key, (Q_con, b_con, c0_con)) in quadratic_lowerthan
        if SparseArrays.nnz(Q_con) != 3 || norm(b_con) > 1e-6 || abs(c0_con) > 1e-6
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
        # TODO continue
    end
end

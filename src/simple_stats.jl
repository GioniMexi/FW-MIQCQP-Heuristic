using LinearAlgebra
using Dates
using PrettyTables

# Calculate the remaining time limit
function remaining_time_limit(start_time::DateTime, current_time_limit::Float64)
    elapsed_time = Dates.value(Dates.now() - start_time) / 1000
    remaining_time = current_time_limit - elapsed_time
    return remaining_time
end

# Measure elapsed time given a start time
function elapsed_time(start_time::DateTime)::Float64
    elapsed_time = Dates.value(Dates.now() - start_time) / 1000
    return elapsed_time
end

# return the number of variables
function num_variables(model)
    return MOI.get(model, MOI.NumberOfVariables())
end

# return the objective sense
function objective_sense(model)
    return MOI.get(model, MOI.ObjectiveSense())
end

# return the objective function type
function objective_function_type(model)
    return MOI.get(model, MOI.ObjectiveFunctionType())
end

"""
Get list of binary and integer variables.
"""
function get_binary_variables(o)
    return MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}())
end

function get_binary_variables_indices(o)
    return getproperty.(get_binary_variables(o), :value)
end

function get_integer_variables(o)
    return MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.Integer}())
end

function get_integer_variables_indices(o)
    return vcat(getproperty.(get_integer_variables(o), :value), getproperty.(get_binary_variables(o), :value))
end

function get_all_variables(o)
    return MOI.get(o, MOI.ListOfVariableIndices())
end

function get_all_variables_indices(o)
    return getproperty.(get_all_variables(o), :value)
end

function problem_statistics(problem_data, heur_params, verbose)
    verbose && println("Start of Statistics\n")

    verbose && println("Problem type: $(problem_data.type)")
    verbose && println()
    # sets of variable indices from the original model.
    num_vars = num_variables(problem_data.optimizer)
    all_vars = Set(1:num_vars)
    binary_set = Set(get_binary_variables_indices(problem_data.optimizer))
    integer_set = Set(get_integer_variables_indices(problem_data.optimizer))

    # continuous are those not in integer_set (binaries are a subset of integer_set)
    continuous_set = setdiff(all_vars, integer_set)

    verbose && println("Number of variables: $num_vars")
    verbose && println("  Binary: $(length(binary_set))")
    verbose && println("  Integer: $(length(integer_set) - length(binary_set))")
    verbose && println("  Continuous: $(length(continuous_set))")
    problem_data.heur_stats.num_vars = num_vars
    problem_data.heur_stats.num_bin_vars = length(binary_set)
    problem_data.heur_stats.num_int_vars = length(integer_set) - length(binary_set)
    problem_data.heur_stats.num_cont_vars = length(continuous_set)

    @assert num_vars == length(integer_set) + length(continuous_set)

    num_quad = 0
    num_bilinear = 0
    num_convex = 0
    num_nonconvex = 0
    # Given a constraint data tuple (Q, b, c0), return the list of involved variable indices.
    function constraint_variable_indices(constraint_data; tol = 1e-12)
        Q, _, _ = constraint_data
        indices = Set{Int}()
        # Process Q (if Q is sparse/dense, ensure proper looping)
        for (i, j, v) in zip(findnz(Q)...)
            if abs(v) > tol
                push!(indices, i)
                push!(indices, j)
            end
        end
        return indices
    end

    # Classify variable types for indices in a constraint.
    function classify_variable_types(idx_set::Set{Int})
        types = String[]
        #count them
        num_binaries = length(idx_set ∩ binary_set)
        num_integers = length(idx_set ∩ integer_set)
        num_continuous = length(idx_set ∩ continuous_set)
        for idx in idx_set
            if idx in binary_set
                push!(types, "binary($num_binaries)")
            elseif idx in integer_set
                push!(types, "integer($num_integers)")
            elseif idx in continuous_set
                push!(types, "continuous($num_continuous)")
            else
                @assert false
            end
        end
        return unique(types)
    end

    # Helper: Analyze one constraint
    function analyze_constraint(kind::String, data, c_idx)
        Q, _, c0 = data

        Q_sym = Symmetric(Q)

        # check if Q is diagonal
        diagQ = isdiag(Q_sym)

        smallest_eig = NaN
        largest_eig = NaN
        conv = NaN
        # Compute the largest eigenvalue or catch error
        try
            eig_res_largest = eigs(Q_sym; nev = 1, which = :LR)
            largest_eig = eig_res_largest[1][1]
        catch
            largest_eig = NaN
        end

        # Compute the smallest eigenvalue or catch error
        try
            eig_res_smallest = eigs(Q_sym; nev = 1, which = :SR)
            smallest_eig = eig_res_smallest[1][1]
        catch
            smallest_eig = NaN
        end

        if isnan(smallest_eig) || isnan(largest_eig)
            conv = false
        else
            conv = smallest_eig >= -1e-6
        end

        idx_set = constraint_variable_indices(data)
        var_types = classify_variable_types(idx_set)

        if kind == "Objective"
            if verbose
                println("  Variable types: ", join(var_types, ","))
                println("  Diagonal terms only: ", diagQ)
                println("  Convex objective: ", conv)
                println("  Objective smallest eigenvalue: ", round(smallest_eig; digits = 2))
                println("  Objective largest eigenvalue: ", round(largest_eig; digits = 2))
            end
        else
            if diagQ
                num_quad += 1
            else
                num_bilinear += 1
            end
            if conv
                num_convex += 1
            else
                num_nonconvex += 1
            end
            if verbose
                println("Constraint[$kind]#$c_idx:")
                println("  Variable types: ", join(var_types, ","))
                println("  Constant:", c0)
                println("  Diagonal terms only: ", diagQ)
                println("  Convex constraint: ", conv)
                println("  Constraint smallest eigenvalue: ", round(smallest_eig; digits = 2))
                println("  Constraint largest eigenvalue: ", round(largest_eig; digits = 2))
            end
        end
    end

    quad_obj = !isempty(problem_data.quad_objective_data[1])

    if verbose
        println("Objective:")
        println("  Objective sense: $(problem_data.sense)")
        println("  Best known solution: $(problem_data.sol_data.best_known_val)")
    end

    # Check if the objective is quadratic or linear
    if quad_obj
        verbose && println("  Objective type: Quadratic")
        if heur_params.verbose_problem_structure_statistics
            analyze_constraint("Objective", problem_data.quad_objective_data, 0)
        end
    else
        verbose && println("  Objective type: Linear")
        @assert length(problem_data.quad_objective_data[2]) == num_vars "$(length(problem_data.quad_objective_data[2])), $(num_vars)"
    end

    if verbose
        println("Quadratic constraints:")
        println("  Number (<=) quadratics: $(length(problem_data.quadratic_lowerthan))")
        println("  Number (==) quadratics: $(length(problem_data.quadratic_equalto))")
    end
    problem_data.heur_stats.num_quad_cons = length(problem_data.quadratic_lowerthan) + length(problem_data.quadratic_equalto)
    problem_data.heur_stats.num_quad_eq_cons = length(problem_data.quadratic_equalto)
    problem_data.heur_stats.num_quad_leq_cons = length(problem_data.quadratic_lowerthan)

    if heur_params.verbose_problem_structure_statistics
        # Process LessEqual constraints
        cidx = 1
        for (ci, data) in problem_data.quadratic_lowerthan
            analyze_constraint("LessEqual", data, cidx)
            cidx += 1
        end

        # Process Equal constraints
        cidx = 1
        for (ci, data) in problem_data.quadratic_equalto
            analyze_constraint("Equal", data, cidx)
            cidx += 1
        end
        if verbose
            println("  Quadratic constraints without bilinear terms: ", num_quad)
            println("  Quadratic constraints with bilinear terms: ", num_bilinear)
            println("  Number of convex constraints: ", num_convex)
            println("  Number of non-convex constraints: ", num_nonconvex)
        end

        problem_data.heur_stats.num_quad_non_bilinear = num_quad
        problem_data.heur_stats.num_quad_bilinear = num_bilinear
        problem_data.heur_stats.num_quad_convex = num_convex
        problem_data.heur_stats.num_quad_non_convex = num_nonconvex
    end

    if problem_data.type == QUBO_BIPARTITE
        ind1, ind2 = problem_data.bipartition_indices
        verbose && println("\nBipartite structure: $(length(ind1)), $(length(ind2))")
    end
    verbose && println("\nMin cover length: $(length(problem_data.cover_indices))")
    verbose && println("Only integers in the min cover: $(is_cover_set_only_integer(problem_data.cover_indices, problem_data.int_vars))\n")

    verbose && println("\nEnd of Statistics")
end

function print_stats(problem_data)

    # Adjust objective values for maximization problems
    best_known_val = problem_data.sol_data.best_known_val
    best_found_val = problem_data.sol_data.best_found_val
    if problem_data.sense == MOI.MAX_SENSE
        best_known_val *= -1
        best_found_val *= -1
    end

    if problem_data.sol_data.num_sols == 0
        relative_gap = 100.0
    elseif iszero(best_found_val) && iszero(best_known_val)
        relative_gap = 0.0
    elseif iszero(best_found_val) || iszero(best_known_val) || (best_found_val * best_known_val < 0)
        # If either is zero (but not both) or their product is negative, they have opposite signs
        relative_gap = 100.0
    else
        denom = max(abs(best_found_val), abs(best_known_val))
        gap = round(abs(best_found_val - best_known_val) / denom * 100; digits = 2)
        relative_gap = min(gap, 100.0)
    end

    combined_matrix = Any[
        "Variable Types" "#" "-"
        "Total Variables" problem_data.heur_stats.num_vars "-"
        "Binary Variables" problem_data.heur_stats.num_bin_vars "-"
        "Integer Variables" problem_data.heur_stats.num_int_vars "-"
        "Continuous Variables" problem_data.heur_stats.num_cont_vars "-"
        "-" "-" "-"
        "Quadratic Constraint Types" "#" "-"
        "All" problem_data.heur_stats.num_quad_cons "-"
        "(==) Constraints" problem_data.heur_stats.num_quad_eq_cons "-"
        "(<=) Constraints" problem_data.heur_stats.num_quad_leq_cons "-"
        "-" "-" "-"
        "Heuristic" "Solutions Found" "Times Improved"
        "RENS" problem_data.heur_stats.num_sols_rens problem_data.heur_stats.num_sols_improved_rens
        "ARENS" problem_data.heur_stats.num_sols_arens problem_data.heur_stats.num_sols_improved_arens
        "RINS" problem_data.heur_stats.num_sols_rins problem_data.heur_stats.num_sols_improved_rins
        "Undercover" problem_data.heur_stats.num_sols_undercover problem_data.heur_stats.num_sols_improved_undercover
        "Alternating" problem_data.heur_stats.num_sols_alternating problem_data.heur_stats.num_sols_improved_alternating
        "-" "-" "-"
        "Statistics" "Value" "-"
        "Best Known" best_known_val "-"
        "Best Found" best_found_val "-"
        "Relative Gap(%)" relative_gap "-"
        "Total Solutions" problem_data.sol_data.num_sols "-"
        "Total Time" elapsed_time(problem_data.time_ref) "-"
    ]
    # Define headers and subheaders
    headers = ["###", "###", "###"]

    # Print the combined table
    pretty_table(combined_matrix; header = headers)
end

"""
Add the solution to the solution data and write the result.
"""
function add_solution(problem_data, curr_obj_val, x)
    problem_data.sol_data.num_sols += 1
    problem_data.sol_data.best_found_val = curr_obj_val
    problem_data.sol_data.best_found_sol = x
    println("FOUND SOLUTION (time, value, thread): $(elapsed_time(problem_data.time_ref)), $(round(curr_obj_val; digits = 4)), $(Threads.threadid())")
    # write_result_solution(problem_data, curr_obj_val)
    write_solution_file(problem_data, curr_obj_val, x)
end

"""
Returns the objective value for x and a boolean indicating if the solution is improving the best found solution.
"""
function evaluate_objective(problem_data::ProblemData, x::Vector{Float64})::Tuple{Float64,Bool}
    Q_obj, c_obj, c0_obj = problem_data.quad_objective_data
    if isempty(Q_obj)
        curr_obj_val = dot(c_obj, x) + c0_obj
    else
        curr_obj_val = 0.5 * dot(x, Q_obj, x) + dot(c_obj, x) + c0_obj
    end

    improving = true
    # Check if the solution improves the best found solution
    if problem_data.sense == MOI.MAX_SENSE
        curr_obj_val = -curr_obj_val
        if curr_obj_val <= problem_data.sol_data.best_found_val
            improving = false
        end
    elseif curr_obj_val >= problem_data.sol_data.best_found_val
        improving = false
    end
    return curr_obj_val, improving
end

"""
Check every solution found by Boscia.
"""
function build_post_heuristic_callback(problem_data::ProblemData, heur_params::HeurParams)
    return function post_heuristics_callback(tree, node, x, origin)
        time_ref = problem_data.time_ref

        curr_obj_val, improving = evaluate_objective(problem_data, x)

        # No need to check the feasibility if the solution is not improving
        if !improving
            return true, elapsed_time(time_ref), tree.root.problem.f(x), x
        end

        if !Boscia.check_linear_feasibility(problem_data.optimizer, x) || !Boscia.is_integer_feasible(tree, x)
            return true, elapsed_time(time_ref), tree.root.problem.f(x), x
        end
        # Check lower-than constraints
        for (key, val) in problem_data.quadratic_lowerthan
            Q_con, b_con, c0_con = val
            constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
            if constraint_value > heur_params.eps_cons
                if heur_params.verbose_violation
                    println("Infeasible for quadratic_lowerthan with violation: $constraint_value")
                end
                return true, elapsed_time(time_ref), tree.root.problem.f(x), x
            end
        end
        # Check equal-to constraints
        for (key, val) in problem_data.quadratic_equalto
            Q_con, b_con, c0_con = val
            constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
            if !isapprox(constraint_value, 0.0; atol = heur_params.eps_cons)
                if heur_params.verbose_violation
                    println("Infeasible for quadratic_equalto with violation: $(abs(constraint_value))")
                end
                return true, elapsed_time(time_ref), tree.root.problem.f(x), x
            end
        end
        add_solution(problem_data, curr_obj_val, x)
        return true, elapsed_time(time_ref), tree.root.problem.f(x), x
    end
end

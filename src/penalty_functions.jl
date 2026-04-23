"""
Return modified objective function with power penalty, to be used in Frank-Wolfe.
"""
function objective_power_penalty(x, problem_data, mu_barrier, power_penalty)
    obj_penalty = 0.0
    for (key, val) in problem_data.quadratic_lowerthan
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        if constraint_value > 0
            obj_penalty += mu_barrier * constraint_value^power_penalty
        end
    end
    for (key, val) in problem_data.quadratic_equalto
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        if abs(constraint_value) > 1e-6
            obj_penalty += mu_barrier * abs(constraint_value)^power_penalty
        end
    end
    return obj_penalty
end

"""
Return modified gradient with power penalty, to be used in Frank-Wolfe.
"""
function gradient_power_penalty!(storage, x, problem_data, mu_barrier, power_penalty)
    for (key, val) in problem_data.quadratic_lowerthan
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        if constraint_value > 0
            tmp = mu_barrier * power_penalty * constraint_value^(power_penalty - 1)
            mul!(storage, Q_con, x, tmp, 1)
            mul!(storage, I, b_con, tmp, 1)
            # grad_constraint = Q_con * x + b_con
            # storage .+= mu_barrier * grad_constraint * power_penalty * constraint_value^(power_penalty - 1)
        end
    end
    for (key, val) in problem_data.quadratic_equalto
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        if abs(constraint_value) > 1e-6
            tmp = mu_barrier * power_penalty * abs(constraint_value)^(power_penalty - 1) * sign(constraint_value)
            mul!(storage, Q_con, x, tmp, 1)
            mul!(storage, I, b_con, tmp, 1)
            # grad_constraint = Q_con * x + b_con
            # storage .+= mu_barrier * grad_constraint * power_penalty * abs(constraint_value)^(power_penalty-1) * sign(constraint_value)
        end
    end
    if isnan(norm(storage))
        @info "Error"
    end
    return nothing
end

"""
Return modified objective function with log barrier penalty, to be used in Frank-Wolfe.
"""
function objective_log_barrier_penalty(x, problem_data, mu_barrier, shift_barrier)
    obj_penalty = 0.0
    # log barrier for leq quadratic constraints
    for (key, val) in problem_data.quadratic_lowerthan
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        obj_penalty -= mu_barrier * log(-constraint_value + shift_barrier)
    end
    return obj_penalty
end

"""
Return modified gradient with log barrier penalty, to be used in Frank-Wolfe.
"""
function gradient_log_barrier_penalty!(storage, x, problem_data, mu_barrier, shift_barrier)
    # Gradient of the logarithm of leq quadratic constraints
    for (key, val) in problem_data.quadratic_lowerthan
        Q_con, b_con, c0_con = val
        constraint_value = 0.5 * dot(x, Q_con, x) + dot(b_con, x) + c0_con
        tmp = mu_barrier / (constraint_value + shift_barrier)
        mul!(storage, Q_con, x, tmp, 1)
        mul!(storage, I, b_con, tmp, 1)
        # grad_constraint = Q_con * x + b_con
        # storage .-= mu_barrier * grad_constraint / (constraint_value + shift_barrier)
    end
    return nothing
end

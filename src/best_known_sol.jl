# A dictionary with "name" => best_sol_value
best_known_solutions = Dict(
)

function get_best_known_solution(prob_name::String)
    # if it is not in the dictionary retun NaN
    return get(best_known_solutions, prob_name, NaN)
end

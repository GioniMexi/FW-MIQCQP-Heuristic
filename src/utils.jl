"""
Return the header of the result file.
"""
function write_result_header(result_path, walltime_start)
    if isempty(result_path) # if result_path is an empty string, do nothing
        return
    end
    if isdir(result_path) # if result_path does exist, empty it from previous solutions
        for entry in readdir(result_path; join = true)
            if occursin("SolutionFile_", entry)
                rm(entry; recursive = true, force = true)
            end
        end
    else # if result_path does not exist, create it
        mkpath(result_path)
    end

    result_file = joinpath(result_path, "ResultFile")

    open(result_file, "w") do io
        println(io, "Walltime_start: $(walltime_start)")
        println(io, "elapsed_walltime_in_seconds\tbest_objective_value")
    end
end

"""
Write a file with the solution.
"""
function write_solution_file(problem_data, new_objective_value::Float64, x::Vector{Float64})
    result_path, result_file, time_ref = problem_data.sol_data.result_path, problem_data.sol_data.result_file, problem_data.time_ref

    # if result_path is an empty string, do nothing
    if isempty(result_path)
        return
    end

    full_result_file = joinpath(result_path, "SolutionFile_t$(Threads.threadid())_s$(problem_data.sol_data.num_sols)")

    # get variable names
    vars = MOI.get(problem_data.optimizer, MOI.ListOfVariableIndices())
    @assert length(MOI.get(problem_data.optimizer, MOI.ListOfVariableIndices())) == length(x)

    x = problem_data.postsolve_function(x)

    # Open the file in write mode and write the new result
    open(full_result_file, "w") do io
        println(io, "$(elapsed_time(time_ref))\t$new_objective_value")
        println(io, "Variable_name \tVariable_value")
        for i in 1:length(x)
            println(io, "$(MOI.get(problem_data.optimizer, MOI.VariableName(), vars[i])) \t $(x[i])")
        end
    end
end

function benchmark_oracles(f, grad!, x_gen, lmo; k = 100, nocache = true)
    x = x_gen()
    sv = sizeof(x) / 1024^2
    println("\nSize of single atom ($(eltype(x))): $sv MB\n")
    to = TimerOutput()
    @showprogress 1 "Testing f... " for i in 1:k
        x = x_gen()
        @timeit to "f" temp = f(x)
    end
    @showprogress 1 "Testing grad... " for i in 1:k
        x = x_gen()
        temp = similar(x)
        @timeit to "grad" grad!(temp, x)
    end
    @showprogress 1 "Testing lmo... " for i in 1:k
        x = x_gen()
        @timeit to "lmo" temp = compute_extreme_point(lmo, x)
    end
    @showprogress 1 "Testing dual gap... " for i in 1:k
        x = x_gen()
        gradient = collect(x)
        grad!(gradient, x)
        v = compute_extreme_point(lmo, gradient)
        @timeit to "dual gap" begin
            dual_gap = FrankWolfe.fast_dot(x, gradient) - FrankWolfe.fast_dot(v, gradient)
        end
    end
    #= @showprogress 1 "Testing update... (Emphasis: OutplaceEmphasis) " for i in 1:k
         x = x_gen()
         gradient = collect(x)
         grad!(gradient, x)
         v = compute_extreme_point(lmo, gradient)
         gamma = 1 / 2
         @timeit to "update (OutplaceEmphasis)" @memory_mode(
             OutplaceEmphasis(),
             x = (1 - gamma) * x + gamma * v
         )
     end
     @showprogress 1 "Testing update... (Emphasis: InplaceEmphasis) " for i in 1:k
         x = x_gen()
         gradient = collect(x)
         grad!(gradient, x)
         v = compute_extreme_point(lmo, gradient)
         gamma = 1 / 2
         # TODO: to be updated to broadcast version once data structure ScaledHotVector allows for it
         @timeit to "update (InplaceEmphasis)" @memory_mode(
             InplaceEmphasis(),
             x = (1 - gamma) * x + gamma * v
         )
     end =#
    if !nocache
        @showprogress 1 "Testing caching 100 points... " for i in 1:k
            @timeit to "caching 100 points" begin
                cache = [gen_x() for _ in 1:100]
                x = gen_x()
                gradient = collect(x)
                grad!(gradient, x)
                v = compute_extreme_point(lmo, gradient)
                gamma = 1 / 2
                test = (x -> FrankWolfe.fast_dot(x, gradient)).(cache)
                v = cache[argmin(test)]
                val = v in cache
            end
        end
    end
    print_timer(to)
    #@show TimerOutputs.time(to["lmo"]), TimerOutputs.time(to["lmo"])/1e9
    #@show TimerOutputs.tottime(to), TimerOutputs.tottime(to)/1e9
    return to
end

"""
Estimate the mu_barrier by computing the 1-norm of the objective function coefficients.
"""
function mu_estimate(quad_objective_data)
    Q, b, c0 = quad_objective_data
    μ = norm(Q, 1) + norm(b, 1)
    return μ
end

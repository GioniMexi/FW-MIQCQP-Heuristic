"""
Watch the new files in a separate thread to write the result file.
"""
function file_watcher(result_path::String, is_minimisation, stop_signal::Threads.Atomic{Bool}, ready_signal::Threads.Event)
    println("File watcher started...")
    seen_files = Set(readdir(result_path)) # Track already seen files
    notify(ready_signal) # Signal that the watcher has started
    result_file = joinpath(result_path, "ResultFile")
    sol_cnt = 1
    best_obj = is_minimisation ? Inf : -Inf

    while !stop_signal[]
        sleep(0.2) # Avoid busy-waiting
        current_files = Set(readdir(result_path))
        new_files = setdiff(current_files, seen_files)
        sleep(0.5) # only activate every half second to avoid taking up CPU

        # Filter only solution files and sort them by the first float in their first line
        solution_files = filter(file -> occursin("SolutionFile_", file), new_files)
        sorted_files = sort(collect(solution_files), by = file -> begin
            full_path = joinpath(result_path, file)
            open(full_path, "r") do f
                first_line = readline(f)
                # Extract the first float from the first line (assumes two floats separated by space)
                return parse(Float64, split(first_line, '\t')[1])
            end
        end)

        for file in sorted_files
            if occursin("SolutionFile_", file)
                full_path = joinpath(result_path, file)
                solution_file = joinpath(result_path, "SolutionFile$sol_cnt")
                open(full_path, "r") do infile
                    open(result_file, "a") do f1
                        first_line = readline(infile)
                        time, obj = parse.(Float64, split(first_line))
                        if (is_minimisation && obj < best_obj) || (!is_minimisation && obj > best_obj)
                            best_obj = obj
                            println("New file detected: $file -> SolutionFile$sol_cnt")
                            println(f1, first_line)
                            # Write remaining lines to second file
                            open(solution_file, "w") do f2
                                for line in eachline(infile)
                                    println(f2, line)
                                end
                            end
                            sol_cnt += 1
                        end
                    end
                end
                # Delete the file after processing
                rm(full_path, force=true)
            end
        end

        seen_files = current_files # Update the seen files
    end

    println("File watcher stopped.")
end

function nthreads_for_type()
    nthreads = Threads.nthreads()
    if nthreads > 1
        return nthreads - 1
    else
        return 1
    end
end

"""
Initialise the parameters depending on the type of the problem.
"""
function initialise_strategy!(h::HeurParams, type::QType, μ)
    h.mu_barrier = μ
    if is_type_qbo(type) # QBO (including QUBO)
        h.node_limit = 1000
    elseif is_type_lmo(type) # LMO (including LBO)
        h.node_limit = 10
        h.max_iter_follow_gradient_heu = 10
    else # QMO
        h.node_limit = 100
        h.max_iter_follow_gradient_heu = 10
    end
end

"""
Modify in place the parameters for the thread depending on the type of the problem and the thread number.
"""
function parallelisation_strategy!(h::HeurParams, type::QType, spectrum::Vector{Float64}, i::Int, nthreads::Int)
    if nthreads > 2 # no strategy for exactly two active threads (so -t 3 with the watchdog) to avoid errors with LinRange
        if i == 1 # first thread always goes for fw_max_iter = 1 (to reduce the time to first solution)
            h.node_limit = 1000
            h.fw_max_iter = 1
        else
            if is_type_qbo(type) # QBO (including QUBO)
                p_min = sum(spectrum .> 0) / length(spectrum)
                h.convexify_objective = LinRange(p_min, 1.0, nthreads - 1)[i-1]
            elseif is_type_lmo(type) # LMO (including LBO)
                h.power_penalty = LinRange(1.2, 1.8, nthreads - 1)[i-1]
                if i != nthreads
                    h.fw_max_iter = rand(h.fw_min_iter:2*h.fw_max_iter)
                end
                if isodd(i)
                    # on even threads, we use the default value
                    h.use_mccormick = false
                end
            else # QMO
                h.fw_max_iter = rand(h.fw_min_iter:2*h.fw_max_iter)
                h.power_penalty = 1.2 + 0.6rand()
            end
        end
    end
end

"""
Modify in place the parameters between restart runs.
"""
function restart_strategy!(h::HeurParams, problem_data::ProblemData, i::Int, nthreads::Int, cnt::Int)
    type = problem_data.type

    # Update mu_barrier based on the adaptive_mu parameter
    if h.adaptive_mu
        # Use solution value to guide mu_barrier if a solution exists
        if isfinite(problem_data.sol_data.best_found_val)
            α = abs(problem_data.sol_data.best_found_val)
            if α ≥ 1e-5
                h.mu_barrier = α
            end
        end
    else
        # Use a simple increasing strategy
        h.mu_barrier += 100.0
        if h.mu_barrier > 10000.0
            h.mu_barrier = 100.0
        end
    end

    if h.num_threads > 1
        if i == 1 # first thread does some random stuff
            h.fw_max_iter = iseven(cnt) + 1 # alternate between 1 and 2 fw_max_iter because why not
            h.power_penalty += 0.05
            if h.power_penalty > 1.8
                h.power_penalty = 1.2
            end
        else
            if is_type_qbo(type) # QBO (including QUBO)
                # probably nothing to change here
            elseif is_type_lmo(type) # LMO (including LBO)
                h.power_penalty = 1.2 + 0.6rand()
            else # QMO
                h.power_penalty = 1.2 + 0.6rand()
            end
        end
    end
end

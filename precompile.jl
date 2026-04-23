using PrimalHeuristic

if !isdir("tmp")
    mkdir("tmp")
end

# to precompile ArgParse
if isempty(ARGS)
    append!(ARGS, ["--file_path", "instances/data/Test_prob_1.mps", "--result_path", "tmp", "--time_limit", "60.0", "--use_mccormick", "true"])
end
PrimalHeuristic.julia_main()

# go through all instances
for i in [1, 2, 3, 6, 8, 11, 14, 15, 27, 36, 42, 54, 92]
    println("Instance $i")
    try
        PrimalHeuristic.main(;
                             file_path = "instances/data/Test_prob_$i.mps",
                             result_path = "tmp/",
                             heur_params = PrimalHeuristic.HeurParams(; time_limit = 30.0)
                            )
    catch e
        println(e)
    end
end

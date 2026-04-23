function parse_command_line()
    arg_parse = ArgParse.ArgParseSettings()

    ArgParse.@add_arg_table! arg_parse begin
        "--file_path"
        help = "path to problem file"
        arg_type = String
        required = true

        "--result_path"
        help = "path to result folder. If empty string then no result file/solutions are written"
        arg_type = String
        required = false
        default = ""

        "--solver_choice"
        help = "solver choice"
        arg_type = String
        required = false
        default = "Gurobi"

        #---------------------------------------------------------------------------
        # Parameter File
        #---------------------------------------------------------------------------
        "--parameter_file"
        help = "The JSON file containing various parameter values."
        arg_type = String
        required = false
        default = ""

        #---------------------------------------------------------------------------
        # Competition Rules
        #---------------------------------------------------------------------------

        "--time_limit"
        help = "time limit for the solver"
        arg_type = Float64
        required = false
        default = 300.0

        "--eps_cons"
        help = "constraint tolerance"
        arg_type = Float64
        required = false
        default = 1e-6

        "--eps_int"
        help = "integer tolerance"
        arg_type = Float64
        required = false
        default = 1e-5

        #---------------------------------------------------------------------------
        # Propagation
        #---------------------------------------------------------------------------

        "--big_m"
        help = "BigM"
        arg_type = Float64
        required = false
        default = 1000000.0

        "--big_m_mccormick"
        help = "BigM McCormick"
        arg_type = Float64
        required = false
        default = 10000.0

        #---------------------------------------------------------------------------
        # LNS Heuristics
        #---------------------------------------------------------------------------

        "--prob_rens"
        help = "Activate RENS heuristic"
        arg_type = Float64
        required = false
        default = 0.0

        "--prob_arens"
        help = "Activate ARENS heuristic"
        arg_type = Float64
        required = false
        default = 1.0

        "--prob_rins"
        help = "Activate RINS heuristic"
        arg_type = Float64
        required = false
        default = 1.0

        "--prob_undercover"
        help = "Activate undercover heuristic"
        arg_type = Float64
        required = false
        default = 1.0

        "--num_covers"
        help = "Number of covers to generate for undercover heuristic"
        arg_type = Int
        required = false
        default = 1

        "--rens_rins_perc_fixed_vars"
        help = "percentage of fixed variables to apply RENS and RINS"
        arg_type = Float64
        required = false
        default = 0.5

        "--use_mccormick"
        help = "activate mccormick relaxation"
        arg_type = Bool
        required = false
        default = false

        "--use_propagation_global"
        help = "activate global propagation"
        arg_type = Bool
        required = false
        default = true

        "--use_propagation_undercover"
        help = "activate undercover propagation"
        arg_type = Bool
        required = false
        default = false

        "--use_external_propagation"
        help = "use external C++ propagation"
        arg_type = Bool
        required = false
        default = true

        #---------------------------------------------------------------------------
        # Boscia
        #---------------------------------------------------------------------------
        "--prob_custom_heu"
        help = "use heuristics associated to the custom Boscia lmo"
        arg_type = Float64
        required = false
        default = 0.0

        "--prob_follow_gradient_heu"
        help = "Activate the follow gradient heuristics"
        arg_type = Float64
        required = false
        default = 1.0

        "--max_iter_follow_gradient_heu"
        help = "Maximum number of iterations in the following gradient heuristic"
        arg_type = Int
        required = false
        default = 100

        "--prob_rounding_01_heu"
        help = "Activate special rounding heuristic for 0/1 polytopes"
        arg_type = Float64
        required = false
        default = 1.0

        "--prob_probability_rounding_heu"
        help = "Activate probability rounding heuristics for 0/1 polytopes"
        arg_type = Float64
        required = false
        default = 1.0

        "--node_limit"
        help = "node limit for Boscia"
        arg_type = Int
        required = false
        default = 100

        "--use_secant"
        help = "Use the Secant line search"
        arg_type = Bool
        required = false
        default = true

        "--fw_lazy"
        help = "Activate lazification in Frank-Wolfe"
        arg_type = Bool
        required = false
        default = true

        "--fw_max_iter"
        help = "Maximum iteration of Frank-Wolfe"
        required = false
        default = 500

        "--fw_min_iter"
        help = "Minimum iteration of Frank-Wolfe"
        required = false
        default = 50

        "--fw_variant"
        help = "Used Frank-Wolfe variant in Boscia (0: default, 1: vanilla)"
        required = false
        default = 0

        "--fw_timeout"
        help = "Safety timeout for Frank-Wolfe"
        arg_type = Float64
        required = false
        default = 15.0

        #---------------------------------------------------------------------------
        # Quadratics Handling
        #---------------------------------------------------------------------------

        "--shift_barrier"
        help = "shift barrier"
        arg_type = Int
        required = false
        default = 1

        "--mu_barrier"
        help = "mu barrier, if zero then we estimate dynamically"
        arg_type = Float64
        required = false
        default = 0.0

        "--adaptive_mu"
        help = "adapt mu dynamically based on solution values"
        arg_type = Bool
        required = false
        default = false

        "--convexify_objective"
        help = "convexify the objective function"
        arg_type = Float64
        required = false
        default = 0.0

        "--power_penalty"
        help = "power penalty for quadratic constraints, log barrier if < 0"
        arg_type = Float64
        required = false
        default = 1.2

        #---------------------------------------------------------------------------
        # Verbosity
        #---------------------------------------------------------------------------

        "--verbose_statistics"
        help = "verbose heuristics"
        arg_type = Bool
        required = false
        default = false

        "--verbose_boscia"
        help = "verbose Boscia"
        arg_type = Bool
        required = false
        default = false

        "--verbose_fw"
        help = "verbose FW"
        arg_type = Bool
        required = false
        default = false

        "--verbose_mip_oracle"
        help = "verbose MIP Oracle"
        arg_type = Bool
        required = false
        default = false

        "--verbose_mip_lns"
        help = "verbose MIP LNS"
        arg_type = Bool
        required = false
        default = false

        "--verbose_violation"
        help = "verbose violation of constraints"
        arg_type = Bool
        required = false
        default = false

        "--verbose_lns"
        help = "verbose statistics of LNS heuristics"
        arg_type = Bool
        required = false
        default = false

        #-----------------------------------------------------------------------------------
        # Restart
        #-----------------------------------------------------------------------------------
        "--restart"
        help = "Restart Boscia"
        arg_type = Bool
        required = false
        default = true

        "--restart_min_time"
        help = "Minimum time in seconds for which doing restart to make sense"
        arg_type = Int
        required = false
        default = 10

        #-------------------------------------------------------------------------------------
        # Fine tuning
        #------------------------------------------------------------------------------------
        "--perform_benchmark_oracle"
        help = "Perform benchmark oracles to estimate LMO cost"
        arg_type = Bool
        required = false
        default = false

        "--prob_alternating"
        help = "Activate alternating heuristic"
        arg_type = Float64
        required = false
        default = 0.0

        "--deperspective"
        help = "Use perspective reformulation when possible"
        arg_type = Bool
        required = false
        default = true

        "--num_threads"
        help = "Number of threads to use (1 thread is for the file watcher)"
        arg_type = Int
        required = false
        default = 7

        "--max_time_lmo"
        help = "Maximum time allowed for each LMO call in seconds"
        arg_type = Float64
        required = false
        default = 10.0
    end

    return ArgParse.parse_args(arg_parse)
end

function julia_main()::Cint
    time_ref = Dates.now()

    @info "Starting parse"
    parsed_args = parse_command_line()

    if isfile(parsed_args["parameter_file"])
        settings_from_file = JSON.Parser.parsefile(parsed_args["parameter_file"])
        for (key, value) in settings_from_file
            parsed_args[key] = value
        end
    end

    heur_params = HeurParams()
    heur_params.solver_choice = parsed_args["solver_choice"]
    heur_params.time_limit = parsed_args["time_limit"]
    heur_params.node_limit = parsed_args["node_limit"]
    heur_params.eps_cons = parsed_args["eps_cons"]
    heur_params.eps_int = parsed_args["eps_int"]
    heur_params.prob_rens = parsed_args["prob_rens"]
    heur_params.prob_arens = parsed_args["prob_arens"]
    heur_params.prob_rins = parsed_args["prob_rins"]
    heur_params.prob_undercover = parsed_args["prob_undercover"]
    heur_params.num_covers = parsed_args["num_covers"]
    heur_params.prob_custom_heu = parsed_args["prob_custom_heu"]
    heur_params.prob_follow_gradient_heu = parsed_args["prob_follow_gradient_heu"]
    heur_params.max_iter_follow_gradient_heu = parsed_args["max_iter_follow_gradient_heu"]
    heur_params.prob_rounding_01_heu = parsed_args["prob_rounding_01_heu"]
    heur_params.prob_probability_rounding_heu = parsed_args["prob_probability_rounding_heu"]
    heur_params.use_secant = parsed_args["use_secant"]
    heur_params.fw_lazy = parsed_args["fw_lazy"]
    heur_params.fw_max_iter = parsed_args["fw_max_iter"]
    heur_params.fw_min_iter = parsed_args["fw_min_iter"]
    heur_params.fw_variant = parsed_args["fw_variant"]
    heur_params.fw_timeout = parsed_args["fw_timeout"]
    heur_params.rens_rins_perc_fixed_vars = parsed_args["rens_rins_perc_fixed_vars"]
    heur_params.use_mccormick = parsed_args["use_mccormick"]
    heur_params.use_propagation_global = parsed_args["use_propagation_global"]
    heur_params.use_propagation_undercover = parsed_args["use_propagation_undercover"]
    heur_params.use_external_propagation = parsed_args["use_external_propagation"]
    heur_params.shift_barrier = parsed_args["shift_barrier"]
    heur_params.mu_barrier = parsed_args["mu_barrier"]
    heur_params.adaptive_mu = parsed_args["adaptive_mu"]
    heur_params.verbose_statistics = parsed_args["verbose_statistics"]
    heur_params.verbose_boscia = parsed_args["verbose_boscia"]
    heur_params.verbose_fw = parsed_args["verbose_fw"]
    heur_params.verbose_mip_oracle = parsed_args["verbose_mip_oracle"]
    heur_params.verbose_mip_lns = parsed_args["verbose_mip_lns"]
    heur_params.verbose_violation = parsed_args["verbose_violation"]
    heur_params.verbose_lns = parsed_args["verbose_lns"]
    heur_params.convexify_objective = parsed_args["convexify_objective"]
    heur_params.power_penalty = parsed_args["power_penalty"]
    heur_params.restart = parsed_args["restart"]
    heur_params.restart_min_time = parsed_args["restart_min_time"]
    heur_params.prob_alternating = parsed_args["prob_alternating"]
    heur_params.deperspective = parsed_args["deperspective"]
    heur_params.perform_benchmark_oracle = parsed_args["perform_benchmark_oracle"]
    heur_params.num_threads = min(7, parsed_args["num_threads"])
    heur_params.max_time_lmo = parsed_args["max_time_lmo"]
    heur_params.bigM = parsed_args["big_m"]
    heur_params.bigM_McCormick = parsed_args["big_m_mccormick"]

    println()
    println(heur_params)

    @info "End parse table"

    main(; time_ref, file_path = parsed_args["file_path"], result_path = parsed_args["result_path"], heur_params)

    return 0
end

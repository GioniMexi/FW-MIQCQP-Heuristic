A Frank-Wolfe-based primal heuristic for Mixed-Integer Quadratically Constrained Programming (MIQCQP).

This code won first place in the [Land-Doig MIP Computational Competition 2025](https://www.mixedinteger.org/2025/).

## Requirements

- **Julia 1.11+** (required)
- **Gurobi 12** with a valid license (required)
- **FixAndProp** (C++ component, optional - see below)

## Compilation

### FixAndProp (C++) - Optional

FixAndProp provides an optimized C++ constraint propagation implementation. If not compiled, global propagation is automatically disabled (`use_propagation_global=false`).

To compile FixAndProp:

```bash
cd FixAndProp
mkdir build && cd build
cmake ..
make -j
```

The compiled library will be at `FixAndProp/build/libfp.{so,dylib,dll}` depending on your platform.

### Julia Setup (requires Julia 1.11+)

1. Unset `GUROBI_HOME` to use the bundled Gurobi version (recommended for reproducibility):

```bash
unset GUROBI_HOME
```

2. Start Julia and activate the project:

```bash
julia --project=.
```

3. Enter package mode (press `]`) and instantiate dependencies:

```
(PrimalHeuristic) pkg> instantiate
```

3. Add the local Boscia and FrankWolfe repository:

```
(PrimalHeuristic) pkg> develop --local ./Boscia.jl
(PrimalHeuristic) pkg> develop --local ./FrankWolfe.jl
```

## Usage

### From Julia REPL

Start Julia with the project:

```bash
julia --project=.
```

Load the module and solve an instance:

```julia
using PrimalHeuristic

# Solve with custom settings and save results
PrimalHeuristic.main(
    file_path="path/to/instance.lp",
    result_path="results/",
    heur_params=PrimalHeuristic.HeurParams(time_limit=30.0, num_threads=7)
)
```

If `result_path=""` then no result files are written.

### Creating a Binary

To create a standalone binary for deployment:

```bash
julia --project=.
```

```julia
using PackageCompiler
create_app(".", "heuristic_binary"; precompile_execution_file="precompile.jl")
```

The binary will be created in the `heuristic_binary/` directory.

### Running with the Binary

```bash
./run_instance.sh <binary_path> <result_dir> <instance> <parameter_file>
```

Or directly:

```bash
./heuristic_binary/bin/PrimalHeuristic \
    --file_path="path/to/instance.mps" \
    --result_path="./results" \
    --parameter_file="params.json"
```

### Parameter File Format

Parameters are specified in JSON format. Example:

```json
{
    "time_limit": 300.0,
    "solver_choice": "Gurobi",
    "num_threads": 1,
    "max_time_lmo": 10.0,
    "fw_max_iter": 100,
    "fw_min_iter": 10,
    "power_penalty": 1.5,
    "mu_barrier": 1000.0,
    "prob_arens": 1.0,
    "prob_rins": 1.0,
    "prob_undercover": 1.0,
    "prob_alternating": 0.0,
    "prob_follow_gradient_heu": 0.0,
    "verbose_boscia": false,
    "verbose_lns": false,
    "verbose_statistics": true
}
```

## Instances

The QPLIB instances used in the paper are available at [QPLIB](https://qplib.zib.de/).

## Citing

If you use this code, please cite:

```bibtex
@misc{mexi2025frankwolfebasedprimalheuristicquadratic,
      title={A Frank-Wolfe-based primal heuristic for quadratic mixed-integer optimization},
      author={Gioni Mexi and Deborah Hendrych and S\'ebastien Designolle and Mathieu Besan\c{c}on and Sebastian Pokutta},
      year={2025},
      eprint={2508.01299},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2508.01299},
}
```

# Decomposition.jl

This package provides a flexible framework for applying (Benders) decomposition algorithms to LP models formulated in `JuMP.jl`, focusing on the application to energy system models with a temporal structure. It includes an annotator for `Calliope` models, but can be extended to other models as well. This repository accompanies the paper _A guide on Benders decomposition in energy system modeling_.

To use this, install the package:

```repl
] add https://github.com/sstroemer/Decomposition.jl
```

## Basic Usage

First, import the necessary packages and create a `JuMP` model:

```julia
import Gurobi
using JuMP: JuMP
using Decomposition: Decomposition
const D = Decomposition
const GRB_ENV = Gurobi.Env()

jump_model = D.jump_model_from_file("models/national_scale_744.mps")
```

> Creating the model is (here) done by reading from an MPS file - you can find this, and other example models, here:  
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17859408.svg)](https://doi.org/10.5281/zenodo.17859408)

Next, we create a `Benders.DecomposedModel`, attaching the `Calliope()` annotator:

```julia
model = D.Benders.DecomposedModel(;
    jump_model = jump_model,
    annotator = D.Calliope(),
    f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
    f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
)
```

We then set various attributes to configure the decomposition algorithm:

```julia
D.set_attribute.(
    model,
    [
        # Set number of total timesteps and requested number of temporal blocks (splits).
        D.Benders.Config.TotalTimesteps(744),
        D.Benders.Config.NumberOfTemporalBlocks(12),

        # Specify that we want to solve the main-model using an IPM.        
        D.Solver.AlgorithmIPM(; model = :main),

        # Activate (multi-cut) optimality cuts, relaxing all sub-models with a penalty.
        D.Benders.OptimalityCutTypeMulti(),
        D.Benders.Sub.RelaxationLinked(; penalty = 1e6),

        # Use the default objective for the main model, and make sure all variables are bounded.
        D.Benders.Main.ObjectiveDefault(),
        D.Benders.Main.VirtualSoftBounds(0.0, 1e6),

        # Set a termination criterion, switch to "info" logging, and disable direct mode.
        D.Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 1000),
        D.Benders.Config.ModelVerbosity(3),
        D.Benders.Config.ModelDirectMode(; enable = false),
    ],
)
```

These attributes are only collected until now, so that the package can automatically decide the best order of execution. Therefore, we need to apply everything now:

```julia
D.finalize!(model)
```

The first time we run `finalize!()` for a certain configuration, it will print a log like:

```text
┌ Info: Begin annotating model
└   annotator = Decomposition.Calliope()
┌ Info: Model annotation successful
└   n_sub_models = 45
┌ Info: Saving annotations to file
└   annotation_file = "national_scale_744_T744_n12.djl.json"
```

This indicates that the annotator has successfully identified 45 sub-models in the original model, based on the temporal decomposition specified. These annotations are also saved to a JSON file, to speed up future runs. The next time you execute the same temporal split on the same model, it will load the annotations from this file instead of re-annotating the model, logging:

```text
┌ Info: Loading annotations from file
└   annotation_file = "national_scale_744_T744_n12.djl.json"
```

If you want to modify certain attributes manually, you can do so before actually running the algorithm. For example, to change some solver parameters of the main-model we can run:

```julia
JuMP.set_attribute(D.Benders.main(model), "BarConvTol", 1e-2)
JuMP.set_attribute(D.Benders.main(model), "OptimalityTol", 1e-2)
JuMP.set_attribute(D.Benders.main(model), "FeasibilityTol", 1e-2)
```

Finally, run the decomposition algorithm:

```julia
while !D.iterate!(model)
end
```

`iterate!` will return `true` once the termination criteria are met - until then, this loop allows you to monitor progress or modify certain parameters on-the-fly. While iterating, the package will report the progress in the terminal, e.g.:

```text
╭──────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ >> Decomposition.jl <<                                              [version::0.1.0  //  algorithm::benders] │
├────────┬───────────────────────────┬───────────────────────────┬───────────────────────────┬─────────────────┤
│        │      objective bound      │      current best gap     │    est. execution time    │   added cuts    │
├────────┼─────────────┬─────────────┼───────────────────────────┤─────────────┬─────────────┤────────┬────────┤
│   iter │       lower │       upper │    absolute │    relative │    wall [s] │     cpu [s] │  feas. │   opt. │
├────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼────────┼────────┤
│      0 │  -4.500e+05 │   6.581e+12 │   6.581e+12 │   1.000e+00 │        6.06 │        6.06 │      0 │     45 │
│      1 │  -4.500e+05 │   4.965e+12 │   4.965e+12 │   1.000e+00 │        6.11 │        6.11 │      0 │     90 │
│      2 │  -4.500e+05 │   4.806e+12 │   4.806e+12 │   1.000e+00 │        6.28 │        6.28 │      0 │    135 │
│      3 │  -1.565e+05 │   4.669e+12 │   4.669e+12 │   1.000e+00 │        6.50 │        6.50 │      0 │    180 │
│      4 │   6.789e+04 │   4.669e+12 │   4.669e+12 │   1.000e+00 │        6.82 │        6.82 │      0 │    225 │
│      5 │   3.107e+05 │   4.669e+12 │   4.669e+12 │   1.000e+00 │        7.14 │        7.14 │      0 │    270 │
...
│     61 │   1.942e+06 │   2.019e+06 │   7.651e+04 │   3.790e-02 │       15.57 │       15.57 │      0 │   2790 │
│     62 │   1.942e+06 │   1.978e+06 │   3.565e+04 │   1.802e-02 │       15.66 │       15.66 │      0 │   2835 │
│     63 │   1.942e+06 │   1.978e+06 │   3.565e+04 │   1.802e-02 │       15.75 │       15.75 │      0 │   2880 │
│     64 │   1.942e+06 │   1.978e+06 │   3.565e+04 │   1.802e-02 │       15.84 │       15.84 │      0 │   2925 │
│     65 │   1.942e+06 │   1.948e+06 │   5.879e+03 │   3.018e-03 │       15.93 │       15.93 │      0 │   2970 │
┌ Info: Termination criterion reached: [opt_gap_rel]
│   iterations = 66
│   seconds = 18.370660466
│   opt_gap_abs = 5879.459616965149
└   opt_gap_rel = 0.003017745799475283
┌ Info: Model statistics
│   lb = 1.942415727468525e6
│   ub = 1.9482951870854902e6
│   ncuts_feas = 0
└   ncuts_opt = 2970
```

After termination, we can summarize the timings collected during the run:

```julia
D.summarize_timings(model)
```

This shows a summary using `TimerOutputs.jl`, e.g.:

```text
Section                ncalls     time    %tot     avg     alloc    %tot      avg
─────────────────────────────────────────────────────────────────────────────────
sub [1 ... 45]          2.97k    8.27s   51.9%  2.79ms    390MiB   54.3%   135KiB
  optimize              2.97k    6.30s   39.5%  2.12ms   60.3MiB    8.4%  20.8KiB
  extract solution      2.97k    1.52s    9.5%   512μs    308MiB   43.0%   106KiB
  fix variables         2.97k    433ms    2.7%   146μs   21.1MiB    2.9%  7.26KiB
main                      132    7.66s   48.1%  58.0ms    328MiB   45.7%  2.48MiB
  optimize                 66    3.16s   19.9%  47.9ms   20.6MiB    2.9%   320KiB
  cuts (generate)          66    3.00s   18.8%  45.4ms    188MiB   26.1%  2.84MiB
  extract solution         66    898ms    5.6%  13.6ms   66.8MiB    9.3%  1.01MiB
  cuts (add)               66    599ms    3.8%  9.07ms   52.7MiB    7.3%   817KiB
  cuts (preprocess)        66    229μs    0.0%  3.48μs   1.03KiB    0.0%    16.0B
  cuts (postprocess)       66    148μs    0.0%  2.24μs   1.03KiB    0.0%    16.0B
─────────────────────────────────────────────────────────────────────────────────
```

## Writing an Annotator

Refer to the `Calliope()` annotator implementation in `src/external/annotators/esm/calliope.jl` for an example of how to write your own annotator for other model types: Simply implementing `_generate_annotations` for a new annotator should be sufficient.

## Citation

<!-- If you use this package in your research, please cite: -->
To be announced.

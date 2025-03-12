# Run once:
# julia --project=experiments experiments/05_bounds.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

const GRB_ENV = Gurobi.Env()
const EXPERIMENT = split(basename(@__FILE__), ".")[1]
const EXPERIMENT_UUID = string(UUIDs.uuid4())
const RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

function experiment(jump_model::JuMP.Model, attributes::Vector; T::Int64, n::Int64)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
        f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
    )

    set_attribute.(
        model,
        [
            Benders.Config.TotalTimesteps(T),
            Benders.Config.NumberOfTemporalBlocks(n),
            Benders.Config.ModelVerbosity(1),
            Benders.Config.ModelDirectMode(; enable = true),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty = 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 1000),
        ],
    )

    set_attribute.(model, attributes)

    finalize!(model)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

attr = Dict(
    :ex1 => [Solver.AlgorithmSimplex(; model = :main)],
    :ex2 => [Solver.AlgorithmSimplex(; model = :main), Benders.Main.VirtualSoftBounds(0.0, 1e6)],
    :ex3 => [Solver.AlgorithmIPM(; model = :main)],
    :ex4 => [Benders.Main.VirtualSoftBounds(0.0, 1e6), Solver.AlgorithmIPM(; model = :main)],
    :ex5 => [Solver.AlgorithmIPM(; model = :main), Benders.Main.RegularizationLevelSet(alpha = 0.5, infeasible_alpha_step = 0.2)],
    :ex6 => [Benders.Main.VirtualSoftBounds(0.0, 1e6), Solver.AlgorithmIPM(; model = :main), Benders.Main.RegularizationLevelSet(alpha = 0.5, infeasible_alpha_step = 0.2)],
)

# Make sure everything's compiled using a small model first.
for (k, v) in attr
    model = experiment(jump_model_from_file("national_scale_120.mps"), v; T = 120, n = 3)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_2184.mps")

# Now run the experiment.
for (k, v) in attr
    println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(k)")
    model = experiment(jump_model, v; T = 2184, n = 12)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "timer_$(k).json"), TimerOutputs.todict(model.timer))
    JSON3.write(
        joinpath(RESULT_DIR, "history_$(k).json"),
        [
            Dict(
                "k" => (entry[:iteration] + 1),
                "t" => entry[:time],
                "lb" => entry[:lower_bound],
                "ub" => entry[:upper_bound],
            ) for entry in model.info[:history]
        ],
    )
end

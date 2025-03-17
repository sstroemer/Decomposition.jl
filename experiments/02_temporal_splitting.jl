# Run N times:
# julia --project=experiments experiments/02_temporal_splitting.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

function experiment(jump_model::JuMP.Model; T::Int64, n::Int64)
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
            Benders.Config.ModelVerbosity(3),
            Benders.Config.ModelDirectMode(; enable = false),
            Solver.AlgorithmIPM(; model = :main),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty = 1e6),
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 250),
        ],
    )

    finalize!(model)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

# Make sure everything's compiled using a small model first.
experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n = 3)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Now run the experiment.
for n in [1, 4, 12, 60, 365]
    println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(n)")
    model = experiment(jump_model; T = 8760, n = n)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "timer_$(n).json"), TimerOutputs.todict(model.timer); allow_inf = true)
    JSON3.write(
        joinpath(RESULT_DIR, "history_$(n).json"),
        [
            Dict(
                "k" => (entry[:iteration] + 1),
                "t" => entry[:time],
                "lb" => entry[:lower_bound],
                "ub" => entry[:upper_bound],
            ) for entry in model.info[:history]
        ];
        allow_inf = true,
    )
end

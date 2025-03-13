# Run N times:
# julia --project=experiments experiments/03_dropping_cuts.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition


GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))


function experiment(jump_model::JuMP.Model; T::Int64, n::Int64, drop::Int64)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
        f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
    )

    if drop > 0
        set_attribute(model, Benders.CutPostprocessingDropNonBinding(; iterations = drop, threshold=1e-8))
    elseif drop == -1
        set_attribute(model, Benders.CutPreprocessingRemoveRedundant(rtol_coeff=1e-3, rtol_const=1e-3))
    elseif drop < -1
        set_attribute(model, Benders.CutPreprocessingRemoveRedundant(rtol_coeff=1e-3, rtol_const=1e-3))
        set_attribute(model, Benders.CutPostprocessingDropNonBinding(; iterations = -drop, threshold=1e-8))
    end

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
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 1000),
        ],
    )

    finalize!(model)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

# Make sure everything's compiled using a small model first.
for drop in [-75, -1, 0, 75]
    experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n = 3, drop)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_2184.mps")

# Now run the experiment.
for drop in [-75, -50, -40, -35, -30, -1, 0]
    model = experiment(jump_model; T = 2184, n = 42, drop)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "timer_$(drop).json"), TimerOutputs.todict(model.timer))
end

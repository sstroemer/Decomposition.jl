# Run N times:
# julia --project=experiments experiments/03_dropping_cuts.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

function experiment(jump_model::JuMP.Model; T::Int64, n::Int64, drop::Int64, ipm::Bool)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
        f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
    )

    # NOTE: 1e-6 is the default FeasibilityTol for Gurobi.

    if drop > 0
        set_attribute(model, Benders.CutPostprocessingDropNonBinding(; iterations = drop, threshold = 1e-6))
    elseif drop == -1
        set_attribute(model, Benders.CutPreprocessingRemoveRedundant(; rtol_coeff = 1e-6, rtol_const = 1e-6))
    elseif drop < -1
        set_attribute(model, Benders.CutPreprocessingRemoveRedundant(; rtol_coeff = 1e-6, rtol_const = 1e-6))
        set_attribute(model, Benders.CutPostprocessingDropNonBinding(; iterations = -drop, threshold = 1e-6))
    end

    ipm && set_attribute(model, Solver.AlgorithmIPM(; model = :main))

    set_attribute.(
        model,
        [
            Benders.Config.TotalTimesteps(T),
            Benders.Config.NumberOfTemporalBlocks(n),
            Benders.Config.ModelVerbosity(3),
            Benders.Config.ModelDirectMode(; enable = false),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty = 1e6),
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 500),
        ],
    )

    finalize!(model)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

# Examples: n=24, drop=[0, -1, 40, 50, 60, 70, 80, 90], ipm=true
n, drop, ipm = length(ARGS) == 3 ? parse.((Int64, Int64, Bool), ARGS) : (60, 50, true)

# Make sure everything's compiled using a small model first.
experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n, drop, ipm)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Run now.
println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(n) / $(drop) / $(ipm)")
model = experiment(jump_model; T = 8760, n, drop, ipm)
JSON3.write(joinpath(RESULT_DIR, "timer_$(n)_$(drop)_$(ipm).json"), TimerOutputs.todict(model.timer); allow_inf = true)

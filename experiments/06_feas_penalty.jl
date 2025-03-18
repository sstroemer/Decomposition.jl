# Run N times:
# julia --project=experiments experiments/06_feas_penalty.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

function experiment(jump_model::JuMP.Model, f; T::Int64, n::Int64, penalty::Float64)
    model = Benders.DecomposedModel(; jump_model, annotator = Calliope(), f_opt_main = f, f_opt_sub = f)

    set_attribute.(
        model,
        [
            Benders.Config.TotalTimesteps(T),
            Benders.Config.NumberOfTemporalBlocks(n),
            Benders.Config.ModelVerbosity(3),
            Benders.Config.ModelDirectMode(; enable = false),  # NOTE: if `true`, line 40 in "relaxation.jl" taking up >> 50% of the overall time.
            Solver.AlgorithmIPM(; model = :main),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty),
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Main.RegularizationLevelSet(; alpha = 0.25, infeasible_alpha_step = 0.25),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 250),
        ],
    )

    finalize!(model)

    # Explicitly set this for the level set.
    if JuMP.solver_name(Benders.main(model)) == "Gurobi"
        JuMP.set_attribute(Benders.main(model), "BarConvTol", 1e-3)
    elseif JuMP.solver_name(Benders.main(model)) == "HiGHS"
        JuMP.set_attribute(Benders.main(model), "ipm_optimality_tolerance", 1e-3)
    end

    while !iterate!(model; nthreads = -1)
    end

    return model
end

# Make sure everything's compiled using a small model first.
experiment(
    jump_model_from_file("national_scale_120.mps"),
    () -> Gurobi.Optimizer(GRB_ENV);
    T = 120,
    n = 3,
    penalty = 1e6,
)
experiment(jump_model_from_file("national_scale_120.mps"), () -> HiGHS.Optimizer(); T = 120, n = 3, penalty = 1e6)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Now run the experiment.
i = length(ARGS) == 1 ? parse(Int64, ARGS[1]) : 7

println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(π)")
model_gurobi = experiment(jump_model, () -> Gurobi.Optimizer(GRB_ENV); T = 8760, n = 24, penalty = 10.0^i)

max_slack = max(
    maximum(maximum(abs.(JuMP.value.(m[:z_neg]))) for m in Benders.subs(model_gurobi)),
    maximum(maximum(abs.(JuMP.value.(m[:z_pos]))) for m in Benders.subs(model_gurobi)),
)

# Cut off if we encounter a slack of at least 1 (MW / ...).
# Setting this too low may otherwise wrongfully trigger for a 0.01 convergence tolerance.
if max_slack >= 1.0
    println("Ecountered non-zero slack: $(max_slack)")
    exit()
end

model_highs = experiment(jump_model, () -> HiGHS.Optimizer(); T = 8760, n = 24, penalty = 10.0^i)

# Write results.
JSON3.write(joinpath(RESULT_DIR, "g_timer_$(i).json"), TimerOutputs.todict(model_gurobi.timer); allow_inf = true)
JSON3.write(joinpath(RESULT_DIR, "h_timer_$(i).json"), TimerOutputs.todict(model_highs.timer); allow_inf = true)

# Run N times:
# julia --project=experiments experiments/06_feas_penalty.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition


const GRB_ENV = Gurobi.Env()
const EXPERIMENT = split(basename(@__FILE__), ".")[1]
const EXPERIMENT_UUID = string(UUIDs.uuid4())
const RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))


function experiment(jump_model::JuMP.Model, f; T::Int64, n::Int64, penalty::Float64)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = f,
        f_opt_sub = f,
    )

    set_attribute.(
        model,
        [
            Benders.Config.TotalTimesteps(T),
            Benders.Config.NumberOfTemporalBlocks(n),
            Benders.Config.ModelVerbosity(1),
            Benders.Config.ModelDirectMode(; enable = false),  # NOTE: if `true`, line 40 in "relaxation.jl" taking up >> 50% of the overall time.
            Solver.AlgorithmIPM(; model = :main),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty),
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Main.RegularizationLevelSet(alpha = 0.1, infeasible_alpha_step = 0.2),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 1000),
        ],
    )

    finalize!(model)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

# Make sure everything's compiled using a small model first.
experiment(jump_model_from_file("national_scale_120.mps"), () -> Gurobi.Optimizer(GRB_ENV); T = 120, n = 3, penalty = 1e6)
# experiment(jump_model_from_file("national_scale_120.mps"), () -> HiGHS.Optimizer(); T = 120, n = 3, penalty = 1e6)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
for i in 0:50
    τ = 0.5
    π = 1e8 * τ^i
    if π < 1.0
        break
    end

    println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(π)")
    model_gurobi = experiment(jump_model, () -> Gurobi.Optimizer(GRB_ENV); T = 744, n = 8, penalty = π)

    max_slack = max(
        maximum(maximum(abs.(JuMP.value.(m[:z_neg]))) for m in Benders.subs(model_gurobi)),
        maximum(maximum(abs.(JuMP.value.(m[:z_pos]))) for m in Benders.subs(model_gurobi)),
    )
    
    if max_slack > 1e-3
        println("Ecountered non-zero slack: $(max_slack)")
        break
    end

    # model_highs = experiment(jump_model, () -> HiGHS.Optimizer(); T = 744, n = 8, penalty = π)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "g_timer_$(i).json"), TimerOutputs.todict(model_gurobi.timer))
    # JSON3.write(joinpath(RESULT_DIR, "h_timer_$(i).json"), TimerOutputs.todict(model_highs.timer))
end

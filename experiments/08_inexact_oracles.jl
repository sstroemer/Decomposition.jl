# Run N times:
# julia --project=experiments experiments/08_inexact_oracles.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

function experiment(jump_model::JuMP.Model; T::Int64, n::Int64, tol::Float64)
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
            Benders.Main.RegularizationLevelSet(; alpha = 0.25, infeasible_alpha_step = 0.25),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 250),
        ],
    )

    finalize!(model)

    # Explicitly set method/crossover/tolerance here to make it verbose.
    cutoff = JuMP.num_variables(model.models[2]) * 0.5
    for sm in Benders.subs(model)
        if (tol > 0) && (JuMP.num_variables(sm) > cutoff)
            JuMP.set_attribute(sm, "Method", 2)
            JuMP.set_attribute(sm, "Crossover", 0)
            JuMP.set_attribute(sm, "BarConvTol", tol)
        else
            JuMP.set_attribute(sm, "Method", 3)
            JuMP.set_attribute(sm, "Crossover", 1)
            JuMP.set_attribute(sm, "BarConvTol", 1e-8)
        end
    end

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
experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n = 2, tol = 0.0)
experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n = 2, tol = 1e-5)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Run now.
i = length(ARGS) == 1 ? parse(Int64, ARGS[1]) : 3
println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(i)")
model = experiment(jump_model; T = 8760, n = 1, tol = (i == 0) ? 0.0 : 10.0^(-i))

# Write results.
JSON3.write(joinpath(RESULT_DIR, "timer_$(i).json"), TimerOutputs.todict(model.timer); allow_inf = true)
JSON3.write(
    joinpath(RESULT_DIR, "history_$(i).json"),
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

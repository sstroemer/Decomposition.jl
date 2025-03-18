# Run once:
# julia --project=experiments experiments/05_bounds.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition

GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))

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
            Benders.Config.ModelVerbosity(3),
            Benders.Config.ModelDirectMode(; enable = false),
            Benders.OptimalityCutTypeMulti(),
            Benders.Sub.RelaxationLinked(; penalty = 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 250),
        ],
    )

    set_attribute.(model, attributes)

    finalize!(model)

    if any(a isa Benders.Main.RegularizationLevelSet for a in attributes)
        # Explicitly set this for the level set.
        if JuMP.solver_name(Benders.main(model)) == "Gurobi"
            JuMP.set_attribute(Benders.main(model), "BarConvTol", 1e-3)
        elseif JuMP.solver_name(Benders.main(model)) == "HiGHS"
            JuMP.set_attribute(Benders.main(model), "ipm_optimality_tolerance", 1e-3)
        end
    end

    while !iterate!(model; nthreads = -1)
    end

    return model
end

attr = Dict(
    :ex1 => [Solver.AlgorithmSimplex(; model = :main)],
    :ex2 => [Solver.AlgorithmSimplex(; model = :main), Benders.Main.VirtualSoftBounds(0.0, 1e6)],
    :ex3 => [Solver.AlgorithmIPM(; model = :main)],
    :ex4 => [Benders.Main.VirtualSoftBounds(0.0, 1e6), Solver.AlgorithmIPM(; model = :main)],
    :ex5 => [
        Solver.AlgorithmIPM(; model = :main),
        Benders.Main.RegularizationLevelSet(; alpha = 0.25, infeasible_alpha_step = 0.25),
    ],
    :ex6 => [
        Benders.Main.VirtualSoftBounds(0.0, 1e6),
        Solver.AlgorithmIPM(; model = :main),
        Benders.Main.RegularizationLevelSet(; alpha = 0.25, infeasible_alpha_step = 0.25),
    ],
)

# Make sure everything's compiled using a small model first.
for (k, v) in attr
    model = experiment(jump_model_from_file("national_scale_120.mps"), v; T = 120, n = 3)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Now run the experiment.
for (k, v) in attr
    println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(k)")
    model = experiment(jump_model, v; T = 8760, n = 24)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "timer_$(k).json"), TimerOutputs.todict(model.timer); allow_inf = true)
    JSON3.write(
        joinpath(RESULT_DIR, "history_$(k).json"),
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

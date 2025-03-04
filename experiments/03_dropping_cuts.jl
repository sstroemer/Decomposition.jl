using JuMP: JuMP
import HiGHS, Gurobi
using Decomposition

const GRB_ENV = Gurobi.Env()


function experiment(jump_model::JuMP.Model, attributes::Vector; T::Int64, n::Int64)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
        f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
    )

    set_attribute.(model, attributes)

    set_attribute.(model, [
        Benders.Config.TotalTimesteps(T),
        Benders.Config.NumberOfTemporalBlocks(n),
        Benders.Config.ModelVerbosity(1),
        Benders.Config.ModelDirectMode(enable=true),
        Solver.AlgorithmIPM(model = :main),
        Benders.OptimalityCutTypeMulti(),
        Benders.Sub.RelaxationLinked(penalty=1e6),
        Benders.Main.VirtualSoftBounds(0.0, 1e6),
        Benders.Main.ObjectiveDefault(),
        Benders.Termination.Stop(opt_gap_rel = 1e-2, iterations = 1000),
    ])

    finalize!(model)

    while !iterate!(model; nthreads = -1); end
    
    return model
end

# Make sure everything's compiled using a small model first.
model = experiment(jump_model_from_file("national_scale_120.mps"), [Benders.CutPostprocessingDropNonBinding(iterations=50)]; T=120, n=3)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
for drop in [30, 90] # 20:10:90
    model = experiment(jump_model, [Benders.CutPostprocessingDropNonBinding(iterations=drop)]; T=744, n=12)

    iter = current_iteration(model)
    time = round(model.timer["main"].accumulated_data.time / 1e6 / iter; digits=2)
    println("drop = $(drop) \t :: \t iterations = $(iter) \t | \t avg_time_main = $(time)ms")
end

# NOTE: `drop <= 15` fails to converge in `<= 500` iterations

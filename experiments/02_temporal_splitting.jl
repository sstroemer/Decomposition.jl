using JuMP: JuMP
import HiGHS, Gurobi
using Decomposition

const GRB_ENV = Gurobi.Env()


function experiment(jump_model::JuMP.Model; T::Int64, n::Int64)
    model = Benders.DecomposedModel(;
        jump_model,
        annotator = Calliope(),
        f_opt_main = () -> Gurobi.Optimizer(GRB_ENV),
        f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV),
    )

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
model = experiment(jump_model_from_file("national_scale_120.mps"); T=120, n=3)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
for n in [8, 12] # [1, 2, 3, 4, 6, 8, 12, 24, 31, 62]
    model = experiment(jump_model; T=744, n=n)

    iter = current_iteration(model)
    time = round(model.timer["main"].accumulated_data.time / 1e9; digits=2)
    println("n = $(n) \t :: \t iterations = $(iter) \t | \t time_main = $(time)s")
end

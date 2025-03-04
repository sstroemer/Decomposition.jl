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
        Benders.Config.ModelDirectMode(enable=false),
        Benders.Config.TotalTimesteps(T),
        Benders.Config.NumberOfTemporalBlocks(n),
        Solver.AlgorithmIPM(model = :main),
        Benders.Main.VirtualSoftBounds(0.0, 1e6),
        Benders.Main.ObjectiveDefault(),
        Benders.Termination.Stop(opt_gap_rel = 1e-2, iterations = 1000),
    ])

    finalize!(model)

    while !iterate!(model; nthreads = -1); end
    
    return model
end

attr = Dict(
    :ex1 => [Benders.OptimalityCutTypeMulti(), Benders.Sub.RelaxationLinked(penalty=1e6)],
    :ex2 => [Benders.OptimalityCutTypeMulti(), Benders.FeasibilityCutTypeMulti(), Solver.ExtractDualRay(model = :sub)],
    :ex3 => [Solver.DualizeModel(model = :sub), Benders.CutTypeMISFSZ()],
)

# Make sure everything's compiled using a small model first.
for attrs in values(attr)
    model = experiment(jump_model_from_file("national_scale_120.mps"), attrs; T=120, n=3)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
model = experiment(jump_model, attr[:ex1]; T=744, n=12)
model = experiment(jump_model, attr[:ex2]; T=744, n=12)
model = experiment(jump_model, attr[:ex3]; T=744, n=12)

summarize_timings(model)

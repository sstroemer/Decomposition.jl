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

        # Solver.AlgorithmIPM(model = :main),
        # Benders.OptimalityCutTypeMulti(),
        # Benders.FeasibilityCutTypeMulti(),
        # Benders.CutTypeMISFSZ(),

        Benders.Main.ObjectiveDefault(),
        # Benders.Main.VirtualSoftBounds(0.0, 1e6),
        Benders.Termination.Stop(opt_gap_rel = 1e-2, iterations = 1000),
    ])

    finalize!(model)

    while !iterate!(model; nthreads = -1); end
    
    return model
end

attr_base = [Benders.OptimalityCutTypeMulti(), Benders.Sub.RelaxationLinked(penalty=1e6)]
attr_spec = Dict(
    :ex1 => [],
    :ex2 => [Solver.AlgorithmIPM(model = :main)],
    :ex3 => [Benders.Main.VirtualSoftBounds(0.0, 1e6)],
    :ex4 => [Solver.AlgorithmIPM(model = :main), Benders.Main.VirtualSoftBounds(0.0, 1e6)],
    # :ex2 => [Benders.OptimalityCutTypeMulti(), Benders.FeasibilityCutTypeMulti(), Solver.ExtractDualRay(model = :sub)],
)

# Make sure everything's compiled using a small model first.
for attrs in values(attr_spec)
    model = experiment(jump_model_from_file("national_scale_120.mps"), [attr_base; attrs]; T=120, n=3)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
model = experiment(jump_model, [attr_base; attr_spec[:ex1]]; T=744, n=12)
model = experiment(jump_model, [attr_base; attr_spec[:ex2]]; T=744, n=12)
model = experiment(jump_model, [attr_base; attr_spec[:ex3]]; T=744, n=12)
model = experiment(jump_model, [attr_base; attr_spec[:ex4]]; T=744, n=12)

summarize_timings(model)

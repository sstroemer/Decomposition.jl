# Run N times:
# julia --project=experiments experiments/07_vi.jl

import TimerOutputs, JSON3, UUIDs
import JuMP, HiGHS, Gurobi
using Decomposition


GRB_ENV = Gurobi.Env()
EXPERIMENT = split(basename(@__FILE__), ".")[1]
EXPERIMENT_UUID = string(UUIDs.uuid4())
RESULT_DIR = mkpath(joinpath(@__DIR__, "out", EXPERIMENT, EXPERIMENT_UUID))


function cb_post_annotate(model; τ::Float64, merge_first::Bool)
    merge_first || τ > 0.0 || return nothing

    n_sub_models = count((k) -> startswith(string(k), "sub_"), keys(model.annotations[:variables]))
    ml = maximum(length(v) for (k, v) in model.annotations[:variables] if startswith(string(k), "sub_"))
    
    model.annotations[:variables][:main_vi_sub] = Int64[]
    model.annotations[:constraints][:main_ci_sub] = Int64[]
    del_sub = []
    for si in 1:n_sub_models
        sm = Symbol("sub_$si")
        if (merge_first && sm == :sub_1) || length(model.annotations[:variables][sm]) < ml * τ
            append!(model.annotations[:variables][:main_vi_sub], model.annotations[:variables][sm])
            append!(model.annotations[:constraints][:main_ci_sub], model.annotations[:constraints][sm])
            push!(del_sub, si)
        end
    end
    
    reduce = 0
    for si in 1:n_sub_models
        if si in del_sub
            sm = Symbol("sub_$si")
            delete!(model.annotations[:variables], sm)
            delete!(model.annotations[:constraints], sm)
            reduce += 1
        elseif reduce > 0
            sm_src = Symbol("sub_$si")
            sm_dst = Symbol("sub_$(si-reduce)")
            model.annotations[:variables][sm_dst] = model.annotations[:variables][sm_src]
            model.annotations[:constraints][sm_dst] = model.annotations[:constraints][sm_src]
            delete!(model.annotations[:variables], sm_src)
            delete!(model.annotations[:constraints], sm_src)
        end
    end

    return nothing
end

function experiment(jump_model::JuMP.Model; T::Int64, n::Int64, τ::Float64, merge_first::Bool, feasibility::Bool = true)
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
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 1000),
        ],
    )

    feasibility && set_attribute(model, Solver.ExtractDualRay(model = :sub))
    feasibility && set_attribute(model, Benders.FeasibilityCutTypeMulti())
    feasibility || set_attribute(model, Benders.Sub.RelaxationLinked(penalty = 1e6))

    finalize!(model; cb_post_annotate = (model) -> cb_post_annotate(model; τ, merge_first))

    while !iterate!(model; nthreads = -1)
    end

    return model
end

attr = Dict(
    :ex1 => Dict(:τ => 0.00, :merge_first => false),
    :ex2 => Dict(:τ => 0.01, :merge_first => false),
    :ex3 => Dict(:τ => 0.10, :merge_first => false),
    :ex4 => Dict(:τ => 0.00, :merge_first => true),
    :ex5 => Dict(:τ => 0.01, :merge_first => true),
    :ex6 => Dict(:τ => 0.10, :merge_first => true),
)

# Make sure everything's compiled using a small model first.
jump_model = jump_model_from_file("national_scale_120.mps")
experiment(jump_model; T = 120, n = 3, τ = 0.0, merge_first = false, feasibility = false)
for (k, v) in attr
    experiment(jump_model; T = 120, n = 3, v...)
end

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_744.mps")

# Now run the experiment.
for (k, v) in attr
    println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(k)")
    model = experiment(jump_model; T = 744, n = 12, v...)

    # Write results.
    JSON3.write(joinpath(RESULT_DIR, "timer_$(k).json"), TimerOutputs.todict(model.timer))
end

println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> baseline")
model = experiment(jump_model; T = 744, n = 12, τ = 0.0, merge_first = false, feasibility = false)
JSON3.write(joinpath(RESULT_DIR, "timer_baseline.json"), TimerOutputs.todict(model.timer))

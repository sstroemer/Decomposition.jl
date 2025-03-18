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

function experiment(jump_model::JuMP.Model; T::Int64, n::Int64, τ::Float64, merge_first::Bool, feasibility::Bool = true, ipm::Bool)
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
            Benders.Main.VirtualSoftBounds(0.0, 1e6),
            Benders.Main.ObjectiveDefault(),
            Benders.Termination.Stop(; opt_gap_rel = 1e-2, iterations = 250),
        ],
    )

    ipm && set_attribute(model, Solver.AlgorithmIPM(; model = :main))

    feasibility && set_attribute(model, Solver.ExtractDualRay(; model = :sub))
    feasibility && set_attribute(model, Benders.FeasibilityCutTypeMulti())
    feasibility || set_attribute(model, Benders.Sub.RelaxationLinked(; penalty = 1e6))

    finalize!(model; cb_post_annotate = (model) -> cb_post_annotate(model; τ, merge_first))

    ipm || JuMP.set_attribute(Benders.main(model), "Method", 3)

    while !iterate!(model; nthreads = -1)
    end

    return model
end

attr = Dict(
    :ex0 => Dict(:τ => 0.00, :merge_first => false, :ipm => true, :feasibility => false),
    :ex1 => Dict(:τ => 0.00, :merge_first => false, :ipm => true),
    :ex2 => Dict(:τ => 0.01, :merge_first => false, :ipm => true),
    :ex3 => Dict(:τ => 0.01, :merge_first => false, :ipm => false),
    :ex4 => Dict(:τ => 0.01, :merge_first => true, :ipm => false),
    
    # :ex3 => Dict(:τ => 0.10, :merge_first => false),
    # :ex4 => Dict(:τ => 0.00, :merge_first => true),
    # :ex5 => Dict(:τ => 0.01, :merge_first => true),
    # :ex6 => Dict(:τ => 0.10, :merge_first => true),
)

exp_idx = length(ARGS) == 1 ? ARGS[1] : "1"
sel = Symbol("ex$(exp_idx)")

# Make sure everything's compiled using a small model first.
experiment(jump_model_from_file("national_scale_120.mps"); T = 120, n = 3, attr[sel]...)

# Load JuMP model.
jump_model = jump_model_from_file("national_scale_8760.mps")

# Now run the experiment.
println("Running experiment: $(EXPERIMENT) >> $(EXPERIMENT_UUID) >> $(exp_idx)")
model = experiment(jump_model; T = 8760, n = 24, attr[sel]...)

# Write results.
JSON3.write(joinpath(RESULT_DIR, "timer_$(exp_idx).json"), TimerOutputs.todict(model.timer); allow_inf = true)

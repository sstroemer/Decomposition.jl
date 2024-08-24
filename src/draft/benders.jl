bd_main(model::DecomposedModel) = model.models[1]
bd_sub(model::DecomposedModel; index::Int = 1) = model.models[1 + index]
bd_subs(model::DecomposedModel) = model.models[2:end]

struct BD_MainObjectiveObj <: DecompositionAttribute end
struct BD_MainObjectiveCon <: DecompositionAttribute end

struct BD_SubObjectivePure <: DecompositionAttribute; index::Int64; end
struct BD_SubObjectiveFull <: DecompositionAttribute; index::Int64; end

struct BD_SubEnsureFeasibility <: DecompositionAttribute
    index::Int64
    penalty::Float64
    penalty_map::Dict{ConstraintRef, AffExpr}

    BD_SubEnsureFeasibility(index::Int64, penalty::Float64) = new(index, penalty, Dict{ConstraintRef, AffExpr}())
    BD_SubEnsureFeasibility(index::Int64) = new(index, 1e8, Dict{ConstraintRef, AffExpr}())
end

function bd_modify(model::DecomposedModel, attribute::DecompositionAttribute)
    @error "Not implemented"
end

function bd_modify(model::DecomposedModel, ::BD_MainObjectiveObj)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]
    
    @objective(m, Min, model.lpmd.c[idx_v]' * m[:x].data + model.lpmd.c_offset)

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectivePure)
    m = bd_sub(model; index=attribute.index)
    idx_v = model.idx_model_vars[1 + attribute.index]

    obj = AffExpr(0.0)
    for i in idx_v
        (i in model.idx_model_vars[1]) && continue
        add_to_expression!(obj, model.lpmd.c[i], m[:x][i])
    end
    @objective(m, Min, obj)

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubEnsureFeasibility)
    m = bd_sub(model; index=attribute.index)

    # Get all variables that are part of the main-model, and find the constraints related to them.
    # idx_v = model.idx_model_vars[1]
    # idx_c_torelax = findall(sum(model.lpmd.A[:, idx_v] .!= 0; dims=2)[:, 1] .!= 0)

    # TODO: only relax the necessary constraints ...
    # relax_with_penalty!(bd_sub(model), Dict(linking_constraints .=> penalty))
    merge!(attribute.penalty_map, relax_with_penalty!(m; default = attribute.penalty))
    
    return nothing
end

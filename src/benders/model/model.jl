include("cache.jl")

Base.Broadcast.broadcastable(model::Benders.DecomposedModel) = Ref(model)

# Dualization is just flagged, and never directly applied as separate step.
apply!(model::Benders.DecomposedModel, attribute::Solver.DualizeModel) = true

# ... same for these.
apply!(model::Benders.DecomposedModel, attribute::Benders.Config.ModelDebug) = true
apply!(model::Benders.DecomposedModel, attribute::Benders.Config.ModelDirectMode) = true

function apply!(model::Benders.DecomposedModel, attribute::Solver.AbstractAttribute)
    models = attribute.model == :main ? [Benders.main(model)] : Benders.subs(model)
    
    ret = true
    for m in models
        ret &= Solver._modify_jump(m, attribute)
    end
    
    return ret
end

function finalize!(model::Benders.DecomposedModel)
    # Construct & apply all annotations if this is called for the first time.
    if isempty(model.models)
        @timeit model.timer "annotations (generate)" annotate!(model)
        @timeit model.timer "annotations (apply)" apply_annotations!(model)
    end

    # Apply all attributes, that are flagged as dirty.
    for ac in model.attributes
        ac.dirty || continue

        ret = @timeit model.timer "apply attributes" apply!(model, ac.attribute)
        if !ret
            @error "Failed to apply attribute" model_type = typeof(model) attribute = ac.attribute
        else
            ac.dirty = false
        end
    end

    return nothing
end

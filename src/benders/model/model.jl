include("cache.jl")
include("conversion/conversion.jl")

main(model::Benders.DecomposedModel) = model.models[1]
sub(model::Benders.DecomposedModel; index::Int) = model.models[1 + index]
subs(model::Benders.DecomposedModel) = model.models[2:end]

Base.Broadcast.broadcastable(model::DecomposedModel) = Ref(model)

function finalize!(model::Benders.DecomposedModel)
    # Construct & apply all annotations if this is called for the first time.
    if isempty(model.models)
        @timeit model.timer "annotations (generate)" annotate!(model)
        @timeit model.timer "annotations (apply)" _apply_annotations!(model)
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

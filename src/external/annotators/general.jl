abstract type AbstractGeneralAnnotator end

function annotate!(model::AbstractDecomposedModel)
    @info "Begin annotating model" annotator = model.annotator
    return _call_annotator(model, model.annotator)
end

function _call_annotator(model::AbstractDecomposedModel, annotator::AbstractGeneralAnnotator)
    @error "Method `annotate!` not implemented for annotator" model annotator
    return false
end

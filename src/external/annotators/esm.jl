abstract type AbstractExternalESMAnnotator <: AbstractGeneralAnnotator end

include("esm/calliope.jl")

function _call_annotator(model::AbstractDecomposedModel, annotator::AbstractExternalESMAnnotator)
    _generate_annotations(model, annotator)

    # Un-dirty all related attributes.
    for ac in model.attributes
        if ac.attribute isa Benders.Config.TotalTimesteps
            ac.dirty = false
        elseif ac.attribute isa Benders.Config.NumberOfTemporalBlocks
            ac.dirty = false
        end
    end

    return true
end

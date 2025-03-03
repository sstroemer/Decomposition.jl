abstract type AbstractExternalESMAnnotator <: AbstractGeneralAnnotator end

include("esm/calliope.jl")

function _call_annotator(model::AbstractDecomposedModel, annotator::AbstractExternalESMAnnotator)
    T = get_attribute(model, Benders.Config.TotalTimesteps, :T)
    n = get_attribute(model, Benders.Config.NumberOfTemporalBlocks, :n)
    annotation_file = "$(model.name)_T$(T)_n$(n).djl.json"

    if isfile(annotation_file)
        @info "Loading annotations from file" annotation_file
        annotations = JSON3.read(annotation_file, Dict{Symbol, Dict{Symbol, Vector{Int64}}})
        merge!(model.annotations, annotations)
    else
        _generate_annotations(model, annotator)
        @info "Saving annotations to file" annotation_file
        JSON3.write(annotation_file, model.annotations)
    end

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

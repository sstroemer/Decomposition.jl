abstract type AbstractExternalESM end

function annotate!(model::AbstractDecomposedModel, ext_fw::AbstractExternalESM)
    @info "Begin annotating model" ESM = typeof(ext_fw)
    @timeit model.timer "annotate" _generate_annotations(model, ext_fw)

    # Un-dirty all related attributes.
    for ac in model.attributes
        if ac.attribute isa Benders.Config.TotalTimesteps
            ac.dirty = false
        elseif ac.attribute isa Benders.Config.NumberOfTemporalBlocks
            ac.dirty = false
        end
    end

    return nothing
end

include("calliope.jl")

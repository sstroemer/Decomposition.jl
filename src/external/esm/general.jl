abstract type AbstractExternalESM end

# TODO: this should be split into "annotate!(...)" and "build!(...)"

function generate_annotation(model::AbstractDecomposedModel, ext_fw::AbstractExternalESM)
    @error "Not implemented"
end

include("calliope.jl")

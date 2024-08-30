abstract type ExternalESM end

function generate_annotation(model::AbstractDecomposedModel, ext_fw::ExternalESM)
    @error "Not implemented"
end

include("calliope.jl")

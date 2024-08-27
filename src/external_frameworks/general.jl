abstract type ExternalFramework end

function generate_annotation(model::DecomposedModel, ext_fw::ExternalFramework)
    @error "Not implemented"
end

include("calliope.jl")

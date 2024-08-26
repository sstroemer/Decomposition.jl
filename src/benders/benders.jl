bd_main(model::DecomposedModel) = model.models[1]
bd_sub(model::DecomposedModel; index::Int = 1) = model.models[1 + index]
bd_subs(model::DecomposedModel) = model.models[2:end]

bd_has_attribute_type(model::DecomposedModel, type::Type{T}) where T <: DecompositionAttribute = any(attr isa type for attr in model.attributes)

function bd_modify(::DecomposedModel, ::DecompositionAttribute)
    @error "Not implemented"
end

include("main.jl")
include("sub.jl")
include("queries.jl")

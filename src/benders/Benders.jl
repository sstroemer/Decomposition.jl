module Benders

import JuMP
import Printf

import ..Decomposition: AbstractDecompositionAttribute, AbstractDecompositionQuery, AbstractDecomposedModel
import ..Solver: AbstractSolverAttribute

include("model.jl")

function has_attribute_type(model::DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    return any(attr isa type for attr in model.attributes)
end

function get_attributes(model::DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    attributes = AbstractDecompositionAttribute[]
    for attr in model.attributes
        attr isa type && push!(attributes, attr)
    end

    if isempty(attributes)
        @error "Model does not contain an attribute of type `$(T)`, check first with `has_attribute_type(...)`"
    end

    return attributes
end

function get_attribute(model::DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    attributes = get_attributes(model, type)

    if length(attributes) > 1
        @error "Attributes of type `$(T)` not unique, use `get_attributes(...)` instead"
        return nothing
    end

    return attributes[1]
end

main(model::DecomposedModel) = model.models[1]
sub(model::DecomposedModel; index::Int) = model.models[1 + index]
subs(model::DecomposedModel) = model.models[2:end]

include("main/Main.jl")
include("sub/Sub.jl")

function modify(model::DecomposedModel, attribute::AbstractSolverAttribute)
    models = attribute.model == :main ? [main(model)] : subs(model)
    for m in models
        modify(m, attribute)
    end
    add_attribute!(model.attributes, attribute)
    return nothing
end

end

import .Benders
const BDModel = Benders.DecomposedModel

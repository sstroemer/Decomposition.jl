module Decomposition

import JuMP
import Statistics: mean
import SparseArrays
import LinearAlgebra
import OrderedCollections: OrderedDict
import Graphs
import JSON3
import SHA
import Printf

abstract type AbstractDecompositionAttribute end
abstract type AbstractDecompositionQuery end
abstract type AbstractDecomposedModel end

function add_attribute!(model::AbstractDecomposedModel, attribute::AbstractDecompositionAttribute)
    push!(model.attributes, attribute)
    push!(model._attribute_iteration_info, current_iteration(model))
    return nothing
end

function modify(model::AbstractDecomposedModel, attribute::AbstractDecompositionAttribute)
    @error "`modify(...)` not implemented for these types" model = typeof(model) attribute = typeof(attribute)
    return nothing
end

function query(model::AbstractDecomposedModel, query::AbstractDecompositionQuery)
    @error "`query(...)` not implemented for these types" model = typeof(model) query = typeof(query)
    return false
end

function has_attribute_type(model::AbstractDecomposedModel, ::Type{T}) where T <: AbstractDecompositionAttribute
    @error "`has_attribute_type(...)` not implemented for these types" model = typeof(model) attribute_type = T
    return false
end

function get_attributes(model::AbstractDecomposedModel, ::Type{T}) where T <: AbstractDecompositionAttribute
    @error "`get_attributes(...)` not implemented for these types" model = typeof(model) attribute_type = T
    return AbstractDecompositionAttribute[]
end

function get_attribute(model::AbstractDecomposedModel, ::Type{T}) where T <: AbstractDecompositionAttribute
    @error "`get_attribute(...)` not implemented for these types" model = typeof(model) attribute_type = T
    return nothing
end

function Base.show(io::IO, attribute::AbstractDecompositionAttribute)
    str = "$(typeof(attribute))("
    str *= join(["$(prop)=$(getfield(attribute, prop))" for prop in propertynames(attribute)], ",")
    str *= ")"
    print(io, str)
    return nothing
end

include("util/util.jl")

include("external/solver/Solver.jl")

include("benders/Benders.jl")

include("external/esm/general.jl")

end

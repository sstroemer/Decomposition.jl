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
import Logging
using TimerOutputs

const MOI = JuMP.MOI

"""
    AbstractDecompositionAttribute

Abstract type for a decomposition attribute, which is a property of a model that is used in the decomposition process.
"""
abstract type AbstractDecompositionAttribute end

"""
    AbstractDecompositionQuery

Abstract type for a decomposition query, which is a property of a model that is used to query results.
"""
abstract type AbstractDecompositionQuery end

"""
    AbstractDecomposedModel

Abstract type for a decomposed model, which is a model that is decomposed into a multiple models.

# Mandatory fields
- `annotator`
- `info`
- `attributes`
- `timer`
- `results`
- `_defaults`
"""
abstract type AbstractDecomposedModel end

"""
    AbstractCut

Abstract type for a cut, which is a constraint that is added to a model to improve the solution quality.
"""
abstract type AbstractCut end

@kwdef mutable struct DecompositionAttributeContainer
    const model::AbstractDecomposedModel
    const attribute::AbstractDecompositionAttribute
    const iteration::Int = current_iteration(model)
    dirty::Bool = true
end

# TODO: move all of this into a MOIWrapper, examples:
# - https://github.com/ds4dm/Tulip.jl/tree/master/src/Interfaces/MOI
# - https://github.com/atoptima/Coluna.jl/blob/master/src/MOIwrapper.jl
# which then also prevents the clash when doing "using JuMP" and "using Decomposition"

function set_attribute(model::AbstractDecomposedModel, attribute::AbstractDecompositionAttribute)
    push!(model.attributes, DecompositionAttributeContainer(; model, attribute))
    return nothing
end

function has_attribute_type(model::AbstractDecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    return any(ac.attribute isa type for ac in model.attributes)
end

function get_attributes(model::AbstractDecomposedModel, type::Type{T}; suppress_error::Bool = true) where T <: AbstractDecompositionAttribute
    attributes = T[]
    for ac in model.attributes
        ac.attribute isa type && push!(attributes, ac.attribute)
    end

    # If no attributes where found, check defaults of the model.
    if isempty(attributes)
        for attr in model._defaults
            if attr isa type
                push!(attributes, attr)
                break
            end
        end
    end

    if isempty(attributes) && !suppress_error
        @error "Model does not contain an attribute of type `$(T)`, check first with `has_attribute_type(...)`"
    end

    return attributes
end

function get_attribute(model::AbstractDecomposedModel, type::Type{T}; suppress_error::Bool = true, get_first::Bool = false) where T <: AbstractDecompositionAttribute
    attributes = get_attributes(model, type; suppress_error)

    isempty(attributes) && return nothing

    if length(attributes) > 1
        @error "Attributes of type `$(T)` not unique, use `get_attributes(...)` instead"
        return nothing
    end

    return attributes[get_first ? 1 : end]::T
end

function get_attribute(model::AbstractDecomposedModel, type::Type{T}, field::Symbol, args...; get_first::Bool = false) where T <: AbstractDecompositionAttribute
    # We need to force `suppress_error = false` here, otherwise a user might think the `nothing` may be the value of the
    # field, instead of the "attribute not found" return value.
    attributes = get_attributes(model, type; suppress_error = false)

    isempty(attributes) && return nothing

    if length(attributes) > 1
        @error "Attributes of type `$(T)` not unique, use `get_attributes(...)` instead"
        return nothing
    end

    attribute = attributes[get_first ? 1 : end]::T
    isempty(args) && return getfield(attribute, field)
    return hasfield(attribute, field) ? getfield(attribute, field) : args[1]
end

function apply!(model::AbstractDecomposedModel, attribute::AbstractDecompositionAttribute)
    @info "No suitable implementation for `apply!`" typeof(model) typeof(attribute)
    return false
end

function execute!(model::AbstractDecomposedModel, query::AbstractDecompositionQuery)
    @info "No suitable implementation for `execute!`" typeof(model) typeof(query)
    return false
end

function cache_get(model::AbstractDecomposedModel, entry::Symbol)
    if !haskey(model._cache, entry)
        model._cache[entry] = getfield(@__MODULE__, Symbol("_cache_build_", entry))(model)    
    end
    
    return model._cache[entry]
end

include("util/util.jl")
include("external/solver/Solver.jl")
include("external/annotators/general.jl")
include("benders/Benders.jl")
include("external/annotators/esm.jl")

import .Benders

# This is based on how JuMP automates exports, with slight modifications.
# See: https://github.com/jump-dev/JuMP.jl which is licensed under MPL v2.0,
# and lists: Copyright (c) 2017: Iain Dunning, Joey Huchette, Miles Lubin, and contributors.
for sym in names(@__MODULE__; all = true)
    (sym in [Symbol(@__MODULE__), :eval, :include]) && continue
    sym_string = string(sym)   
    (startswith(sym_string, "_") || startswith(sym_string, "@_")) && continue
    Base.isidentifier(sym) || (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end])) || continue
    @eval export $sym
end

end

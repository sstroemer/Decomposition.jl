"""
    AbstractBendersDecompositionNode

Mandatory fields: `main`, `sub`
"""
abstract type AbstractBendersDecompositionNode <: AbstractDecompositionNode end

@kwdef struct SimpleBendersDecompositionNode <: AbstractBendersDecompositionNode
    id::ID
    parent::AbstractModelNode
    children::Vector{AbstractNode} = Vector{AbstractNode}()

    main::Vector{AbstractNode} = Vector{AbstractNode}()
    sub::Vector{AbstractNode} = Vector{AbstractNode}()

    ext::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

@kwdef struct DecomposedModelNode <: AbstractModelNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode} = Vector{AbstractNode}()
    model::AbstractModel

    ext::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

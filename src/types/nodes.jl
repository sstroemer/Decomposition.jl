struct Node{T <: AbstractNode}
    Node{T}(parent::AbstractNode; kwargs...) where {T <: AbstractNode} = add_child(parent, T; kwargs...)
end

@kwdef struct NodeRoot <: AbstractNode
    id::UUIDs.UUID
    _uids::Set{UUIDs.UUID} = Set{UUIDs.UUID}()
    children::Vector{AbstractNode} = Vector{AbstractNode}()
end
function Node{NodeRoot}()
    node = NodeRoot(id=UUIDs.uuid4())
    push!(node._uids, node.id)
    return node
end

@kwdef struct NodeModelPresolved <: AbstractModelNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}
    model::AbstractModel

    filename_full::String
    filename_reduced::String
    filename_postsolve::String
end

@kwdef struct NodeFile <: AbstractFileNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}

    filename::String
end

@kwdef struct NodeModel <: AbstractModelNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}
    model::AbstractModel
end

@kwdef struct NodeLP <: AbstractProgramNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}
    program::AbstractProgram
end

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

@kwdef struct FileNodePresolved <: AbstractFileNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}

    filename::String
    filename_original::String
    filename_postsolve::String
end

@kwdef struct FileNodeGeneral <: AbstractFileNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}

    filename::String
end

@kwdef struct ModelNodeGeneral <: AbstractModelNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}
    model::AbstractModel
end

@kwdef struct ProgramNodeLP <: AbstractProgramNode
    id::UUIDs.UUID
    parent::AbstractNode
    children::Vector{AbstractNode}
    program::ProgramLP
end
Base.show(io::IO, node::ProgramNodeLP) = print(io, "$(_to_str(node)): $(node.program.m) × $(node.program.n) (nz = $(length(node.program.A.nzval) + length(node.program.A.nzval)))")

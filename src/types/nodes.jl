struct Node{T <: AbstractNode}
    Node{T}(parent::AbstractNode; kwargs...) where {T <: AbstractNode} = add_child(parent, T; kwargs...)
end

@kwdef struct NodeRoot <: AbstractNode
    id::ID
    _uids::Set{ID} = Set{ID}()
    _optimizer::Any = HiGHS.Optimizer   # TODO: make this configurable
    children::Vector{AbstractNode} = Vector{AbstractNode}()
end
function Node{NodeRoot}()
    node = NodeRoot(id=ID(value = 0))
    push!(node._uids, node.id)
    return node
end
root_node(node::NodeRoot) = node

@kwdef struct FileNodePresolved <: AbstractFileNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode}

    filename::String
    filename_original::String
    filename_postsolve::String

    filename_reduced_sol::String
    filename_full_sol::String
end

@kwdef struct FileNodeGeneral <: AbstractFileNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode}

    filename::String
end

@kwdef struct ModelNodeGeneral <: AbstractModelNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode}
    model::AbstractModel
end

@kwdef struct ModelNodeDualization <: AbstractModelNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode}
    model::AbstractModel

    _primal_model::AbstractModel        # = `parent.model`
    _dual_model::AbstractModel          # = `.model`
    primal_dual_map::Dualization.PrimalDualMap
end

@kwdef struct ProgramNodeLP <: AbstractProgramNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode}
    program::ProgramLP
end
Base.show(io::IO, node::ProgramNodeLP) = print(io, "$(_to_str(node)): $(node.program.m) × $(node.program.n) (nz = $(length(node.program.A.nzval) + length(node.program.A.nzval)))")

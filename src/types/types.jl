# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Internal IDs ├───                                                                                              ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

@kwdef struct ID
    value::Int64
end
Base.show(io::IO, id::ID) = print(io, "$(string(id.value; base=16, pad=5))")
Base.isequal(id1::ID, id2::ID) = id1.value == id2.value
Base.hash(id::ID, h::UInt) = hash(id.value, hash(id.value, h))
Base.isless(id1::ID, id2::ID) = id1.value < id2.value

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Abstract base types: Programs & Nodes ├───                                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

"""
    AbstractNode
    
Mandatory fields: `id`, `_uids`, `parent` (can be `nothing` if node is root), `children`
"""
abstract type AbstractNode end

"""
    AbstractFileNode

Mandatory fields: `filename`
"""
abstract type AbstractFileNode <: AbstractNode end

"""
    AbstractModelNode

Mandatory fields: `model`
"""
abstract type AbstractModelNode <: AbstractNode end

"""
    AbstractProgramNode

Mandatory fields: `program`
"""
abstract type AbstractProgramNode <: AbstractNode end

abstract type AbstractProgram end

abstract type AbstractGeneralModel end
const AbstractModel = Union{AbstractGeneralModel, JuMP.AbstractModel, MOI.ModelLike}

"""
    AbstractLink

Mandatory fields: `from`, `to`, `info
"""
abstract type AbstractLink end

function Base.getproperty(node::AbstractNode, field::Symbol)
    hasproperty(node, field) && return getfield(node, field)
    field == :_uids && return getfield(node, :parent)._uids
    field == :_optimizer && return getfield(node, :parent)._optimizer

    if field in [:parent, :model, :program]
        @debug "$(nameof(typeof(node))).$(field) -> nothing" node = getfield(node, :id)
        return nothing
    end

    @critical "Trying to access unknown field" node = getfield(node, :id) field
end

_to_str(node::AbstractNode) = "$(replace(string(nameof(typeof(node))), "Node" => "")) [$(node.id)]"
Base.show(io::IO, node::AbstractNode) = print(io, _to_str(node))
Base.show(io::IO, node::AbstractFileNode) = print(io, "$(_to_str(node)): $(node.filename)")
Base.show(io::IO, node::AbstractModelNode) = print(io, "$(_to_str(node)): $(JuMP.num_variables(node.model)) vars and $(JuMP.num_constraints(node.model; count_variable_in_set_constraints = false)) cons")

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Base functionality for Nodes ├───                                                                              ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

function _next_id(node::AbstractNode)
    uids = node._uids
    id = ID(maximum(keys(uids)).value + 1)
    return id, uids
end

function add_child(parent::AbstractNode, TypeChildNode::Type{T}; kwargs...) where {T <: AbstractNode}
    id, uids = _next_id(parent)
    child = TypeChildNode(id=id, parent=parent, children=Vector{AbstractNode}(); kwargs...)

    uids[id] = child
    push!(parent.children, child)

    return child
end

root_node(node::AbstractNode) = root_node(node.parent)
node_from_id(node::AbstractNode, id::ID) = get(node._uids, id, nothing)

function solve!(node::AbstractNode)
    @debug "solve!(::AbstractNode)" node
    length(node.children) > 1 && (@critical "Cannot solve multi-children node without additional information" node)

    return (isempty(node.children) ? _solve!(node) : solve!(node.children[1]))
end
_solve!(node::AbstractNode) = @critical "Cannot find implementation of `solve`" node

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Base functionality for Programs ├───                                                                           ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

# TODO

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Detailed types & their functionality ├───                                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

include("programs.jl")
include("nodes.jl")
include("solutions.jl")

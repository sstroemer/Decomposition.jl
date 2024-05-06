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
const AbstractModel = Union{AbstractGeneralModel, JuMP.AbstractModel}

function Base.getproperty(node::AbstractNode, field::Symbol)
    hasproperty(node, field) && return getfield(node, field)
    field == :_uids && return getfield(node, :parent)._uids

    if field in [:parent, :model, :program]
        @debug "$(nameof(typeof(node))).$(field) -> nothing" node = getfield(node, :id)
        return nothing
    end

    @critical "Trying to access unknown field" node = getfield(node, :id) field
end

_to_str(node::AbstractNode) = "$(replace(string(nameof(typeof(node))), "Node" => "")) [$(string(node.id)[1:8])]"
Base.show(io::IO, node::AbstractNode) = print(io, _to_str(node))
Base.show(io::IO, node::AbstractFileNode) = print(io, "$(_to_str(node)): $(node.filename)")
Base.show(io::IO, node::AbstractModelNode) = print(io, "$(_to_str(node)): $(JuMP.num_variables(node.model)) vars, and $(JuMP.num_constraints(node.model; count_variable_in_set_constraints = false)) cons")

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Base functionality for Nodes ├───                                                                              ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

function add_child(parent::AbstractNode, TypeChildNode::Type{T}; kwargs...) where {T <: AbstractNode}
    id = UUIDs.uuid4()

    _uids = parent._uids
    (id in _uids) && @critical "Unexpected ID clash" parent id
    push!(_uids, id)

    child = TypeChildNode(id=id, parent=parent, children=Vector{AbstractNode}(); kwargs...)
    push!(parent.children, child)
    return child
end

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Base functionality for Programs ├───                                                                           ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

# TODO

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ╟───┤ Detailed types & their functionality ├───                                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

include("programs.jl")
include("nodes.jl")

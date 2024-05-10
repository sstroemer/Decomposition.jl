abstract type AbstractResult end 

struct ResultPrimalObjective <: AbstractResult end
struct ResultDualObjective <: AbstractResult end
struct ResultPrimalValue <: AbstractResult end
struct ResultDualValue <: AbstractResult end


const ObjectIndex = Union{MOI.VariableIndex, MOI.ConstraintIndex}
const AbstractObjectIndex = Union{String, MOI.VariableIndex, MOI.ConstraintIndex}
const OptionalObjectIndex = Union{Nothing, MOI.VariableIndex, MOI.ConstraintIndex}
const OptionalAbstractObjectIndex = Union{Nothing, String, MOI.VariableIndex, MOI.ConstraintIndex}


_resolve_object(model::JuMP.Model, index::MOI.VariableIndex) = JuMP.VariableRef(model, index)
_resolve_object(model::JuMP.Model, index::MOI.ConstraintIndex) = JuMP.constraint_ref_with_index(model, index)
_resolve_object(model::JuMP.Model, index::String) = something(JuMP.variable_by_name(model, index), JuMP.constraint_by_name(model, index))


function _get_result_abstract(node::AbstractNode, object::OptionalObjectIndex, type::AbstractResult)
    @debug "get_result(::AbstractNode, ::OptionalAbstractObjectIndex, ::AbstractResult)" node object type

    isempty(node.children) && return _get_result_self(node, object, type)
    length(node.children) > 1 && (@critical "Cannot get result of multi-children node without additional information" node type)
    
    return _get_result_child(node, node.children[1], object, type)
end

get_result(node::AbstractNode, object::OptionalObjectIndex, type::AbstractResult) = _get_result_abstract(node, object, type)
get_result(node::AbstractNode, type::AbstractResult) = get_result(node, nothing, type)
get_result(node::NodeRoot, object::String, type::AbstractResult) = get_result(node.children[1], object, type)

_get_result_child(node::AbstractNode, child::AbstractNode, object::OptionalObjectIndex, type::AbstractResult) = get_result(child, object, type)
_get_result_self(node::AbstractNode, object::OptionalObjectIndex, type::AbstractResult) = @critical "Cannot find implementation of `get_result`" node object type


include("model.jl")
include("file.jl")

@kwdef Solution
    node::AbstractNode

    obj_value_primal::Float64
    obj_value_dual::Float64

    primal_value::Dict{JuMP.VariableRef, Float64}
    reduced_cost::Dict{JuMP.VariableRef, Float64}
    shadow_price::Dict{JuMP.ConstraintRef, Float64}
end

function get_solution(node::AbstractNode)::Solution
    @debug "get_solution(::AbstractNode)" node

    isempty(node.children) && return _get_solution(node)::Solution
    length(node.children) > 1 && (@critical "Cannot get solution of multi-children node without additional information" node)
    
    return get_solution(node.children[1])::Solution
end
_get_solution(node::AbstractNode) = @critical "Cannot find implementation of `get_solution`" node



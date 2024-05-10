function get_result(node::AbstractFileNode, object::String, type::AbstractResult)
    length(node.children) == 1 ||
        @critical "Getting result of file node requries exactly one `FileNodeGeneral` as child" node object type

    child = node.children[1]
    child isa ModelNodeGeneral ||
        @critical "Getting result of file node requries exactly one `FileNodeGeneral` as child" node object type

    jump_object = something(JuMP.variable_by_name(child.model, object), JuMP.constraint_by_name(child.model, object))
    return get_result(child, JuMP.index(jump_object), type)
end

function get_result(node::AbstractFileNode, object::ObjectIndex, type::AbstractResult)
    @warn "Consistency of MOI indices across file writing is currently not ensured" node maxlog = 1
    return _get_result_abstract(node, object, type)
end

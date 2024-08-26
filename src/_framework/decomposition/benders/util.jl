function add_annotation(node::SimpleBendersDecompositionNode; variable::MOI.VariableIndex, annotation::Symbol = :main)
    # TODO: also "remember" for each variable to which stage it belongs to simplify the lookup afterwards

    data = _get_decomposition_data!(node.parent.model)
    annotations = get!(data, :benders_annotations, Dict{Symbol, Set{MOI.VariableIndex}}())
    list = get!(annotations, annotation, Set{MOI.VariableIndex}())
    push!(list, variable)
    return nothing
end

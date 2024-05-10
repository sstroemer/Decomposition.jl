function _solve!(node::AbstractModelNode)
    JuMP.set_optimizer(node.model, node._optimizer)
    JuMP.optimize!(node.model)
    return nothing
end

function solve!(node::ModelNodeDualization)
    @debug "solve!(::ModelNodeDualization)" node
    length(node.children) > 1 && (@critical "Cannot solve multi-children node without additional information" node)

    ret = isempty(node.children) ? _solve!(node) : solve!(node.children[1])
    _update_mock_solution(node)
    return ret
end

function _update_mock_solution(node::ModelNodeDualization)
    # Check if we optimized internally, or need to get results from a child.
    self_opt = isempty(node.children)

    mock_modellike = JuMP.backend(node._outer_model).optimizer.dual_problem.dual_model.model.optimizer
    index_map = node._reference_map_outer_to_model.index_map

    nm_backend = JuMP.backend(node.model)
    mo = node._mock_optimizer

    var_primal(_vi) = self_opt ? MOI.get(nm_backend, MOI.VariablePrimal(), index_map[_vi]) : get_result(node.children[1], index_map[_vi], ResultPrimalValue())
    con_dual(_ci) = self_opt ? MOI.get(nm_backend, MOI.ConstraintDual(), index_map[_ci]) : get_result(node.children[1], index_map[_ci], ResultDualValue())

    primal = [var_primal(vi) for vi in MOI.get(mock_modellike, MOI.ListOfVariableIndices())]
    constraint_duals = [
        ctype => [con_dual(ci) for ci in MOI.get(mock_modellike, MOI.ListOfConstraintIndices{ctype...}())]
        for ctype in MOI.get(mock_modellike, MOI.ListOfConstraintTypesPresent())
    ]
    MOI.Utilities.mock_optimize!(mo, MOI.OPTIMAL, primal, MOI.FEASIBLE_POINT, constraint_duals...)
    
    # Manually set the proper objective values, which allows disabling "internal" calculation of those, which takes a lot of time.
    obj_val = self_opt ? JuMP.objective_value(node.model) : get_result(node.children[1], ResultPrimalObjective())
    dual_obj_val = self_opt ? JuMP.dual_objective_value(node.model) : get_result(node.children[1], ResultDualObjective())
    MOI.set(mo, MOI.ObjectiveValue(), obj_val)
    MOI.set(mo, MOI.DualObjectiveValue(), dual_obj_val)

    return nothing
end

_get_result_self(node::AbstractModelNode, ::Nothing, ::ResultPrimalObjective) = JuMP.objective_value(node.model)
_get_result_self(node::AbstractModelNode, ::Nothing, ::ResultDualObjective) = JuMP.dual_objective_value(node.model)

_get_result_self(node::AbstractModelNode, object::ObjectIndex, ::ResultPrimalValue) = JuMP.value(_resolve_object(node.model, object))
_get_result_self(node::AbstractModelNode, object::ObjectIndex, ::ResultDualValue) = JuMP.dual(_resolve_object(node.model, object))

_get_result_self(node::ModelNodeDualization, ::Nothing, ::ResultPrimalObjective) = JuMP.objective_value(node._outer_model)
_get_result_self(node::ModelNodeDualization, ::Nothing, ::ResultDualObjective) = JuMP.dual_objective_value(node._outer_model)
_get_result_child(node::ModelNodeDualization, child::AbstractNode, ::Nothing, ::ResultPrimalObjective) = JuMP.objective_value(node._outer_model)
_get_result_child(node::ModelNodeDualization, child::AbstractNode, ::Nothing, ::ResultDualObjective) = JuMP.dual_objective_value(node._outer_model)

_get_result_self(node::ModelNodeDualization, object::ObjectIndex, ::ResultPrimalValue) = JuMP.value(_resolve_object(node._outer_model, object))
_get_result_self(node::ModelNodeDualization, object::ObjectIndex, ::ResultDualValue) = JuMP.dual(_resolve_object(node._outer_model, object))
_get_result_child(node::ModelNodeDualization, child::AbstractNode, object::ObjectIndex, ::ResultPrimalValue) = JuMP.value(_resolve_object(node._outer_model, object))
_get_result_child(node::ModelNodeDualization, child::AbstractNode, object::ObjectIndex, ::ResultDualValue) = JuMP.dual(_resolve_object(node._outer_model, object))

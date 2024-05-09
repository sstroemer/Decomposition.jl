function _solve!(node::AbstractModelNode)
    JuMP.set_optimizer(node.model, node._optimizer)
    JuMP.optimize!(node.model)

    return nothing
end

function _get_solution(node::AbstractModelNode)
    if !JuMP.is_solved_and_feasible(node.model)
        @error "Model is not properly solved (unsolved, infeasible, unbounded, or non-optimal)" node
        return nothing
    end

    if !isempty(JuMP.all_constraints(model, JuMP.VariableRef, MOI.Interval{Float64}))
        @critical "Getting solution of ranged constraints is currently not possible" node
    end

    return Solution(
        node=node,
        obj_value_primal=JuMP.objective_value(node.model),
        obj_value_dual=JuMP.dual_objective_value(node.model),
        primal_value=Dict(v => JuMP.value(v) for v in JuMP.all_variables(node.model)),
        reduced_cost=Dict(v => JuMP.reduced_cost(v) for v in JuMP.all_variables(node.model)),
        shadow_price=Dict(c => JuMP.dual(c) for c in JuMP.all_constraints(node.model; include_variable_in_set_constraints=false))
    )
end

function get_solution(node::ModelNodeDualization)::Solution
    length(node.children) > 1 && (@critical "Cannot get solution of multi-children node without additional information" node)
    length(node.children) == 1 && return get_solution(node.children[1])::Solution


    if mode == :primal
        # return Dict(child.primal_dual_map[_getvarname(var)] => value for (var, value) in solution)
    elseif mode == :dual
    end

    @critical "Unknown mode" mode
end

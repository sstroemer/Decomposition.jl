dualize(node::AbstractNode) = dualize(to_model(node))

function dualize(node::AbstractModelNode)
    mock_opt = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}(); eval_objective_value=false, eval_dual_objective_value=false, eval_variable_constraint_dual=false)
    outer_model = JuMP.Model(Dualization.dual_optimizer(() -> mock_opt); add_bridges=false)
    refmap_outer = JuMP.GenericReferenceMap(outer_model, MOI.copy_to(JuMP.backend(outer_model), JuMP.backend(node.model)))

    # Call `optimize!` to create the model inside the mock-optimizer.
    JuMP.optimize!(outer_model)

    # Get the "inner problem" (which is the actual dual problem we want to make accessible for later nodes).
    inner_model = JuMP.Model()
    refmap_inner = JuMP.GenericReferenceMap(outer_model, MOI.copy_to(JuMP.backend(inner_model), JuMP.backend(outer_model).optimizer.dual_problem.dual_model.model.optimizer))
    # inner_model = JuMP.backend(outer_model).optimizer.dual_problem.dual_model

    return Node{ModelNodeDualization}(node, model=inner_model, _mock_optimizer=mock_opt, _reference_map_parent_to_outer=refmap_outer, _outer_model=outer_model, _reference_map_outer_to_model=refmap_inner)

    # dual_problem = Dualization.DualProblem{Float64}(MOI.Utilities.UniversalFallback(DualizableModel{T}()), dual_optimizer)
    # mo = _make_mock_optimizer()
    # mdo = _make_mock_dualoptimizer(mo)
    # dm = JuMP.Model(mdo)
    # dp = Dualization.DualProblem(JuMP.backend(dm))
    # Dualization.dualize(JuMP.backend(node.model), dp)

    # dual_node = Node{ModelNodeDualization}(node, model=dm, _primal_model=node.model, _dual_model=dm, primal_dual_map=dp.primal_dual_map, _mock_optimizer=mo)
    # _update_mock_dualoptimizer(dual_node)
    
    # return dual_node
end

_make_mock_optimizer() = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
_make_mock_dualoptimizer(mo) = Dualization.DualOptimizer{Float64}(MOI.instantiate(MOI.Utilities.Model{Float64}))

function _update_mock_dualoptimizer(node::ModelNodeDualization)
    JuMP.optimize!(node.model)
    # MOI.optimize!(node._mock_optimizer)
end

dualize(node::AbstractNode) = dualize(to_model(node))

function dualize(node::AbstractModelNode)
    dm = JuMP.Model()
    dp = Dualization.DualProblem(JuMP.backend(dm))
    Dualization.dualize(JuMP.backend(node.model), dp)
    return Node{ModelNodeDualization}(node, model=dm, _primal_model=node.model, _dual_model=dm, primal_dual_map=dp.primal_dual_map)
end

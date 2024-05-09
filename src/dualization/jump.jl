dualize(node::AbstractNode) = dualize(to_model(node))

function dualize(node::AbstractModelNode)
    child = Node{ModelNodeGeneral}(node, model=Dualization.dualize(node.model))
    Link{DualLink}(node, child)
    return child
end

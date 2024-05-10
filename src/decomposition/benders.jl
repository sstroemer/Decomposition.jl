"""
    AbstractDecomposition

Mandatory fields: `main`, `sub`
"""
abstract type AbstractBendersDecompositionNode <: AbstractDecompositionNode end

@kwdef struct SimpleBendersDecompositionNode <: AbstractBendersDecompositionNode
    id::ID
    parent::AbstractNode
    children::Vector{AbstractNode} = Vector{AbstractNode}()

    main::Vector{AbstractNode} = Vector{AbstractNode}()
    sub::Vector{AbstractNode} = Vector{AbstractNode}()
end

function SimpleBendersDecomposition(node::AbstractModelNode)
    benders = Node{SimpleBendersDecompositionNode}(node)

    # TODO: Implement this
    push!(benders.main, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))

    # TODO: Implement this
    push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))

    return benders
end

SimpleBendersDecomposition(node::AbstractNode) = SimpleBendersDecomposition(to_model(node))

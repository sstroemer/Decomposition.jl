to_model(node::AbstractNode) = @critical "Failed to convert to model" node

to_model(node::AbstractModelNode) = node

function to_model(node::NodeLP)
    @debug "to_model(::NodeLP)" node
    program = node.program

    # TODO: retain names
    model = JuMP.Model()
    JuMP.@variable(model, x[i = 1:program.n] >= program.lb[i])
    JuMP.@constraint(model, program.A * x .<= program.b)
    JuMP.@constraint(model, program.G * x .== program.h)
    JuMP.@objective(model, Min, SparseArrays.dot(program.c, x) + program.d)

    return Node{NodeModel}(node, model=model)
end

function to_model(node::AbstractFileNode)
    @debug "to_model(::AbstractFileNode)" node
    return Node{NodeModel}(node, model=JuMP.read_from_file(node.filename))
end

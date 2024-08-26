to_model(node::AbstractNode) = @critical "Failed to convert to model" node

to_model(node::AbstractModelNode) = node

function to_model(node::ProgramNodeLP)
    @debug "to_model(::NodeLP)" node
    program = node.program

    # TODO: retain names
    model = JuMP.Model()
    JuMP.@variable(model, x[i = 1:program.n] >= program.lb[i])
    JuMP.@constraint(model, program.A * x .<= program.b)
    JuMP.@constraint(model, program.G * x .== program.h)
    JuMP.@objective(model, Min, SparseArrays.dot(program.c, x) + program.d)

    return Node{ModelNodeGeneral}(node; model = model)
end

function to_model(node::AbstractFileNode)
    @debug "to_model(::AbstractFileNode)" node
    # TODO: remove this after fixing problems in dualize (and "to_program")
    model = JuMP.read_from_file(node.filename)
    _resolve_ranged_constraints!(model)
    return Node{ModelNodeGeneral}(node; model = model)
end

to_model(node::String) = to_model(from_file(node))

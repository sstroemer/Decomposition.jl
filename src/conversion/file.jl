from_file(filename::String) = Node{NodeFile}(Node{NodeRoot}(); filename=filename)

function to_file(node::AbstractNode, filename::String = "")
    @debug "to_file(::AbstractNode, ::String)" node filename
    node = to_model(node)
    (filename == "") && (filename = "$(tempname()).mps")
    JuMP.write_to_file(node.model, filename)
    return Node{NodeFile}(node; filename=filename)
end

function to_file(node::AbstractModelNode, filename::String = "")
    @debug "to_file(::AbstractModelNode, ::String)" node filename
    (filename == "") && (filename = "$(tempname()).mps")
    JuMP.write_to_file(node.model, filename)
    return Node{NodeFile}(node; filename=filename)
end

to_file(node::AbstractFileNode, filename::String = "") = filename == "" ? node : to_file(node, filename)

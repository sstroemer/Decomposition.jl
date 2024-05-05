abstract type AbstractWorkflow end


@kwdef struct Workflow <: AbstractWorkflow
    root::NodeRoot
end

function print_tree(workflow::Workflow; max_depth::Int = -1)
    _walk(workflow.root; max_depth=max_depth)
end

function _walk(node::AbstractNode, depth::Int = 0; visited::Union{Nothing, Set{UUIDs.UUID}} = nothing, max_depth::Int = -1)
    max_depth > 0 && depth > max_depth && return
    visited = something(visited, Set{UUIDs.UUID}())
    node.id in visited && return

    siblings = isnothing(node.parent) ? [] : node.parent.children
    if depth == 0
        println(node)
        println("╤════════════════════════════════════════════")
    elseif !isempty(siblings)
        order = sort([c.id for c in siblings])
        char = node.id == order[end] ? "╰" : "├"
        println("$(repeat(" ", (depth-1) * 3))$(char)─ $(node)")
    else
        println("$(repeat(" ", (depth-1) * 3))$(char)─ $(node)")
    end

    push!(visited, node.id)

    for child in node.children
        _walk(child, depth+1; visited=visited, max_depth=max_depth)
    end
end

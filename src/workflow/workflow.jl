abstract type AbstractWorkflow end


@kwdef struct Workflow <: AbstractWorkflow
    root::NodeRoot
end
Workflow(node::AbstractNode) = Workflow(root=root(node))

function print_tree(workflow::Workflow; max_depth::Int = -1)
    _walk(workflow.root; max_depth=max_depth)
end

function _walk(node::AbstractNode, depth::Int = 0; visited::Union{Nothing, Set{ID}} = nothing, max_depth::Int = -1, open_levels::Union{Nothing, Vector{Int}} = nothing)
    max_depth > 0 && depth > max_depth && return
    visited = something(visited, Set{ID}())
    node.id in visited && return

    siblings = isnothing(node.parent) ? [] : node.parent.children
    has_siblings = length(siblings) > 1
    is_last_child = isempty(siblings) ? false : (node.id == siblings[end].id)

    open_levels = something(open_levels, Vector{Int}())
    is_last_child && !isempty(open_levels) && pop!(open_levels)

    if depth == 0
        println(node)
        println("╤═══════════")
    else
        for i in 1:(depth-1)
            if (i - 1) in open_levels
                print("┆ ")
            else
                print("  ")
            end
        end
        # print("$(repeat(" ", (depth-1) * 2))")
        if has_siblings
            prefix = is_last_child ? "╰" : "├"
            infix = isempty(node.children) ? "─" : "┬"
            println("$(prefix)─$(infix)─ $(node)")
        else
            prefix = "╰"
            infix = isempty(node.children) ? "─" : "┬"
            println("$(prefix)─$(infix)─ $(node)")
        end
    end

    push!(visited, node.id)
    !isempty(node.children) && push!(open_levels, depth)

    for child in node.children
        _walk(child, depth+1; visited=visited, max_depth=max_depth, open_levels=open_levels)
    end
    # pop!(open_levels)
end

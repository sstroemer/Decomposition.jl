include("benders/types.jl")
include("benders/util.jl")


function SimpleBendersDecomposition(node::AbstractModelNode)
    benders = Node{SimpleBendersDecompositionNode}(node)

    # TODO: set optimize hook

    # TODO: Implement this
    # push!(benders.main, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))

    # TODO: Implement this
    # push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    # push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    # push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))
    # push!(benders.sub, Node{ModelNodeGeneral}(benders; model = JuMP.copy_model(node.model)[1]))

    return benders
end

SimpleBendersDecomposition(node::AbstractNode) = SimpleBendersDecomposition(to_model(node))

function finalize(node::SimpleBendersDecompositionNode)
    model_full = node.parent.model
    data = _get_decomposition_data!(model_full)
    annotations = data[:benders_annotations]

    levels = keys(annotations)
    length(levels) in [1, 2] || @critical "SimpleBendersDecomposition currently supports exactly two levels"

    main = annotations[:main]

    model_main, refmap_main = @suppress JuMP.copy_model(model_full)
    vi = first(main)
    con_with_var = collect(c for c in JuMP.all_constraints(model_main; include_variable_in_set_constraints = false) if JuMP.normalized_coefficient(c, JuMP.VariableRef(model_main, vi)) != 0)
    for c in con_with_var
        length(JuMP.constraint_object(c).func.terms) == 1 && continue   # TODO: check for more than 1 var
        JuMP.delete(model_main, c)
    end
    JuMP.delete(model_main, collect(refmap_main[v] for v in JuMP.all_variables(model_full) if JuMP.index(v) ∉ main))
    # TODO
    con_all_main = JuMP.all_constraints(model_main; include_variable_in_set_constraints = false)
    JuMP.delete(model_main, collect(c for c in con_all_main if isempty(JuMP.constraint_object(c).func.terms)))

    node_main = Node{DecomposedModelNode}(node; model=model_main)
    node_main.ext[:refmap] = refmap_main
    push!(node.main, node_main)

    model_sub, refmap_sub = @suppress JuMP.copy_model(model_full)
    node_sub = Node{DecomposedModelNode}(node; model=model_sub)
    node_sub.ext[:refmap] = refmap_sub
    push!(node.sub, node_sub)

    # Create "main -> sub" reference map.
    # refmap_main_to_sub = JuMP.GenericReferenceMap(model_main, ...)
    # TODO
end

function solve!(node::SimpleBendersDecompositionNode)
    finalize(node)

    model_main = node.main[1].model
    model_sub = node.sub[1].model
    var_decomposed = JuMP.all_variables(model_main)
    
    result_lb = -Inf
    result_ub = +Inf
    
    JuMP.@variable(node.main[1].model, θ)
    JuMP.fix(θ, 0.0; force=true)
    JuMP.@objective(node.main[1].model, Min, JuMP.objective_function(node.main[1].model) + θ)
    
    JuMP.set_silent(model_main)
    JuMP.set_silent(model_sub)
    
    println("══════════════════════════════════════════════════")
    println("       [Iterative Benders Decomposition]          ")
    println("─────────┬─────────────┬─────────────┬────────────")
    println("    iter │          lb │          ub │     rel_gap")
    println("─────────┼─────────────┼─────────────┼────────────")
    for k in 1:10
        # Solve main
        solve!(node.main[1])
        obj_value_main = JuMP.objective_value(model_main)
    
        # Get `x_k` and fix it in the subproblem.
        var_value_main = JuMP.value.(var_decomposed)
        for v in var_decomposed
            JuMP.fix(JuMP.VariableRef(model_sub, JuMP.index(v)), JuMP.value(v); force=true)
        end
    
        # Solve sub.
        solve!(node.sub[1])
        obj_value_sub = JuMP.objective_value(model_sub)
    
        # Update bounds, gap, and more.
        result_lb = max(result_lb, obj_value_main)
        result_ub = min(result_ub, (obj_value_main - JuMP.value(θ)) + obj_value_sub)
        gap_abs = result_ub - result_lb
        gap_rel = gap_abs / result_ub
    
        Decomposition._print_iteration(k, result_lb, result_ub, gap_rel)
        if gap_rel < 1e-6 || gap_abs < -Inf
            break
        end
    
        # Get `π_k`.
        π = JuMP.reduced_cost.(collect(JuMP.VariableRef(model_sub, JuMP.index(v)) for v in var_decomposed))
    
        # Add a new cut to the main problem.
        JuMP.is_fixed(θ) && JuMP.unfix(θ)
        JuMP.@constraint(model_main, θ >= obj_value_sub + π' * (var_decomposed .- var_value_main))
    end
    println("══════════════════════════════════════════════════")
end

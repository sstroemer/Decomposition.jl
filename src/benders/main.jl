struct BD_MainObjectiveConstOnly <: DecompositionAttribute end
struct BD_MainObjectiveObj <: DecompositionAttribute end
struct BD_MainObjectiveCon <: DecompositionAttribute end

abstract type BD_AbstractMainCutType <: DecompositionAttribute end
struct BD_MainCutTypeSingle <: BD_AbstractMainCutType end
struct BD_MainCutTypeMulti <: BD_AbstractMainCutType end

struct BD_MainCutType <: BD_AbstractMainCutType
    # options: :off, :single, :multi, :aggregated, :adaptive (changes cuts)

    feasibility::Symbol
    optimality::Symbol
end

struct MainVirtualBounds <: DecompositionAttribute
    lower::Float64
    upper::Float64
end

# TODO: refactor `push!(model.attributes, attribute)` into
# function add_attribute!(model, attribute)
#     push!(model.attributes, attribute) # <=== remember here at WHICH ITERATION it was added
#     return nothing
# end

function bd_modify(model::DecomposedModel, attribute::BD_MainCutType)
    setup_θ = !bd_has_attribute_type(model, BD_MainCutType)

    push!(model.attributes, attribute)

    if setup_θ
        # First time configuring the cut type => setup the θ variable.
        @variable(bd_main(model), θ[i = 1:(length(model.models) - 1)])
        set_lower_bound.(θ, -1e4)           # TODO
    end

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainCutTypeSingle)
    if bd_has_attribute_type(model, BD_MainCutTypeMulti)
        # TODO: allow "changing" the cut type by removing the attribute instead of erroring
        @error "Cannot use both single and multi cuts"
        return nothing
    end
    
    push!(model.attributes, attribute)
    _bd_setup_θ(model)

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainCutTypeMulti)
    if bd_has_attribute_type(model, BD_MainCutTypeSingle)
        # TODO: allow "changing" the cut type by removing the attribute instead of erroring
        @error "Cannot use both single and multi cuts"
        return nothing
    end

    push!(model.attributes, attribute)
    _bd_setup_θ(model)

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::MainVirtualBounds)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]

    for i in idx_v
        isfinite(attribute.lower) && !has_lower_bound(m[:x][i]) && set_lower_bound(m[:x][i], attribute.lower)
        isfinite(attribute.upper) && !has_upper_bound(m[:x][i]) && set_upper_bound(m[:x][i], attribute.upper)
    end

    push!(model.attributes, attribute)
    return nothing
end

function _bd_setup_θ(model::DecomposedModel)
    if !(bd_has_attribute_type(model, BD_MainCutTypeSingle) || bd_has_attribute_type(model, BD_MainCutTypeMulti))
        @error "No valid cut mode (single / multi) specified"
        return nothing
    end

    if haskey(bd_main(model), :θ)
        @error "Variable θ already exists in the main-model"
        return nothing
    end

    if bd_has_attribute_type(model, BD_MainCutTypeSingle)
        @variable(bd_main(model), θ)
        set_lower_bound(θ, -1e4)           # TODO
    else
        @variable(bd_main(model), θ[i = 1:(length(model.models) - 1)])
        set_lower_bound.(θ, -1e4)           # TODO
    end

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveConstOnly)
    m = bd_main(model)
    
    if bd_has_attribute_type(model, BD_MainCutTypeSingle)
        @objective(m, Min, model.lpmd.c_offset + m[:θ])
    else
        @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    end

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveObj)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]

    if bd_has_attribute_type(model, BD_MainCutTypeSingle)
        @objective(m, Min, model.lpmd.c[idx_v]' * m[:x].data + model.lpmd.c_offset + m[:θ])
    else
        @objective(m, Min, model.lpmd.c[idx_v]' * m[:x].data + model.lpmd.c_offset + sum(m[:θ]))
    end

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveCon)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]

    @expression(m, obj, model.lpmd.c[idx_v]' * m[:x].data)

    if bd_has_attribute_type(model, BD_MainCutTypeSingle)
        @objective(m, Min, model.lpmd.c_offset + m[:θ])
    else
        @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    end

    push!(model.attributes, attribute)
    return nothing
end

function bd_generate_cuts(model::DecomposedModel, current_solution::JuMP.Containers.DenseAxisArray)
    # TODO: get the current solution directly from the model

    if !(bd_has_attribute_type(model, BD_MainCutTypeSingle) || bd_has_attribute_type(model, BD_MainCutTypeMulti))
        @error "No valid cut mode (single / multi) specified"
        return nothing
    end

    m_main = bd_main(model)
    m_main_x = m_main[:x]
    vis_main = model.idx_model_vars[1]

    # Track number of added cuts.
    nof_added_opt_cuts = 0
    nof_added_feas_cuts = 0

    # This my or may not be used.
    single_cut_exp = AffExpr(0.0)
    single_opt_cut_empty = true

    new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    for i in 1:(length(model.models) - 1)
        vis_sub = model.idx_model_vars[1 + i]
        m_sub = bd_sub(model; index=i)
        m_sub_x = m_sub[:x]

        cut_type = bd_check_cut_type(m_sub)

        if cut_type == :optimality
            exp_cut = AffExpr(0.0)
            add_to_expression!(exp_cut, model.info[:results][:subs][i][:obj])  # TODO: this should be the "lower bound of the sub-model"

            # TODO: check if this works as expected
            if bd_has_attribute_type(model, BD_MainObjectiveCon)
                add_to_expression!(exp_cut, m_main[:obj])
            end

            for vi in vis_main
                (vi in vis_sub) || continue

                λ = reduced_cost(m_sub_x[vi])
                add_to_expression!(exp_cut, λ, m_main_x[vi])
                add_to_expression!(exp_cut, λ, -current_solution[vi])
            end

            push!(new_cuts[:optimality], (i, exp_cut))
        elseif cut_type == :feasibility
            # Prepare the cut expression.
            exp_cut = AffExpr(model.info[:results][:subs][i][:obj_dual])  # TODO: what's the best way to access this?
            for vi in vis_main
                (vi in vis_sub) || continue

                λ = dual(FixRef(m_sub_x[vi]))
                add_to_expression!(exp_cut, λ, m_main_x[vi])
                add_to_expression!(exp_cut, λ, -current_solution[vi])  # = -fix_value(m_sub_x[vi])
            end

            push!(new_cuts[:feasibility], (i, exp_cut))
        else
            @error "Cut type for sub-model $(i) could not be determined"
        end
    end

    return new_cuts
end

function bd_add_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    return _bd_add_cuts(model, new_cuts, bd_get_attribute(model, BD_AbstractMainCutType))
end

function _bd_add_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, attribute::BD_MainCutTypeSingle)
    nof_added_feas_cuts = 0
    nof_added_opt_cuts = 0

    # TODO: what happens if we aggregate the feasibility cuts too?
    # if !isempty(new_cuts[:feasibility])
    #     exp_agg_cut = AffExpr(0.0)
    #     for (_, exp_cut) in new_cuts[:feasibility]
    #         add_to_expression!(exp_agg_cut, exp_cut)
    #     end

    #     push!(
    #         model.cuts[:feasibility],
    #         @constraint(bd_main(model), exp_agg_cut <= 0)
    #     )
    #     nof_added_feas_cuts += 1
    # end

    for (_, exp_cut) in new_cuts[:feasibility]
        push!(
            model.cuts[:feasibility],
            @constraint(bd_main(model), exp_cut <= 0)
        )
        nof_added_feas_cuts += 1
    end

    # TODO: aggregating cuts is basically "averaging" them; we could: (1) weight the average, (2) cluster them before aggregating to generate `N >= 1` cuts
    if !isempty(new_cuts[:optimality])
        exp_agg_cut = AffExpr(0.0)
        for (_, exp_cut) in new_cuts[:optimality]
            add_to_expression!(exp_agg_cut, exp_cut)
        end

        n = length(new_cuts[:optimality])
        push!(
            model.cuts[:optimality],
            @constraint(bd_main(model), bd_main(model)[:θ] >= exp_agg_cut)
        )
        nof_added_opt_cuts += 1
    end

    return nof_added_feas_cuts, nof_added_opt_cuts
end

function _bd_add_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, attribute::BD_MainCutTypeMulti)
    nof_added_feas_cuts = 0
    nof_added_opt_cuts = 0

    for (_, exp_cut) in new_cuts[:feasibility]
        push!(
            model.cuts[:feasibility],
            @constraint(bd_main(model), exp_cut <= 0)
        )
        nof_added_feas_cuts += 1
    end

    for (i, exp_cut) in new_cuts[:optimality]
        push!(
            model.cuts[:optimality],
            @constraint(bd_main(model), bd_main(model)[:θ][i] >= exp_cut)
        )
        nof_added_opt_cuts += 1
    end

    # TODO: is there a multi- vs. single-cut version of feasibility cuts?

    # bd_has_attribute_type(model, BD_MainCutTypeSingle)
    # bd_has_attribute_type(model, BD_MainCutTypeMulti)
    # cut = @constraint(m_main, m_main[:θ][i] >= exp_cut)
    # push!(model.cuts[:optimality], cut)

    # cut = @constraint(m_main, exp_cut < = 0)  

    # TODO: it can happen that we add a lot of `\theta >= 0` cuts here, due to the "trivial" sub-models.
    #       => catch linear dependent cuts and do not add them.

    # if !single_opt_cut_empty
    #     # Some sub-models added to the single optimality cut, so we need to add it to the main-model.
    #     if bd_has_attribute_type(model, BD_MainObjectiveCon)
    #         single_cut_exp += m_main[:obj]
    #     end

    #     cut = @constraint(m_main, m_main[:θ] >= single_cut_exp)
    #     push!(model.cuts[:optimality], cut)
    #     nof_added_opt_cuts += 1
    # end

    return nof_added_feas_cuts, nof_added_opt_cuts
end

function bd_check_cut_type(model::JuMP.Model; verbose::Bool = true)
    if is_solved_and_feasible(model)
        return :optimality
    end

    if dual_status(model) != MOI.INFEASIBILITY_CERTIFICATE
        verbose && (@error "Turn off presolve, or any setting blocking extraction of dual rays")
        return :error
    end

    return :feasibility
end

struct BD_MainObjectiveConstOnly <: DecompositionAttribute end
struct BD_MainObjectiveObj <: DecompositionAttribute end
struct BD_MainObjectiveCon <: DecompositionAttribute end

struct BD_MainCutTypeSingle <: DecompositionAttribute end
struct BD_MainCutTypeMulti <: DecompositionAttribute end

struct MainVirtualBounds <: DecompositionAttribute
    lower::Float64
    upper::Float64
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

    for i in 1:(length(model.models) - 1)
        vis_sub = model.idx_model_vars[1 + i]
        m_sub = bd_sub(model; index=i)
        m_sub_x = m_sub[:x]

        cut_type = bd_check_cut_type(m_sub)

        if cut_type == :optimality
            exp_cut = bd_has_attribute_type(model, BD_MainCutTypeSingle) ? single_cut_exp : AffExpr(0.0)
            add_to_expression!(exp_cut, model.info[:results][:subs][i][:obj])  # TODO: this should be the "lower bound of the sub-model"

            for vi in vis_main
                (vi in vis_sub) || continue

                λ = reduced_cost(m_sub_x[vi])
                add_to_expression!(exp_cut, λ, m_main_x[vi])
                add_to_expression!(exp_cut, λ, -current_solution[vi])
            end
    
            # Only add the cut here, if it is a multi-cut.
            if bd_has_attribute_type(model, BD_MainCutTypeMulti)
                if bd_has_attribute_type(model, BD_MainObjectiveCon)
                    exp_cut += m_main[:obj]
                end
        
                # TODO: it can happen that we add a lot of `\theta >= 0` cuts here, due to the "trivial" sub-models.
                #       => catch linear dependent cuts and do not add them.
                cut = @constraint(m_main, m_main[:θ][i] >= exp_cut)
                push!(model.cuts[:optimality], cut)
                nof_added_opt_cuts += 1
            else
                single_opt_cut_empty = false
            end
        elseif cut_type == :feasibility
            # TODO: is there a multi- vs. single-cut version of feasibility cuts?

            # Prepare the cut expression.
            exp_cut = AffExpr(0.0)
            for vi in vis_main
                (vi in vis_sub) || continue

                λ = dual(FixRef(m_sub_x[vi]))
                add_to_expression!(exp_cut, λ, fix_value(m_sub_x[vi]))  # TODO: get that from the current solution instead
                add_to_expression!(exp_cut, -λ, m_main_x[vi])
            end

            # Add the new cut to the main-model, and to the list of all cuts.
            cut = @constraint(m_main, exp_cut >= model.info[:results][:subs][i][:obj_dual])  
            push!(model.cuts[:feasibility], cut)
            nof_added_feas_cuts += 1
        else
            @error "Cut type for sub-model $(i) could not be determined"
        end
    end

    if !single_opt_cut_empty
        # Some sub-models added to the single optimality cut, so we need to add it to the main-model.
        if bd_has_attribute_type(model, BD_MainObjectiveCon)
            single_cut_exp += m_main[:obj]
        end

        cut = @constraint(m_main, m_main[:θ] >= single_cut_exp)
        push!(model.cuts[:optimality], cut)
        nof_added_opt_cuts += 1
    end

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

function generate_cuts(model::Benders.DecomposedModel)
    all_new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    _generate_cuts_from_dual(model, all_new_cuts)
    _generate_cuts_from_primal(model, all_new_cuts)

    # Remove redundant (duplicate) cuts.
    # TODO: also, do not create ones that already exist

    new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    for cut_type in [:feasibility, :optimality]
        for cut1 in all_new_cuts[cut_type]
            unique_cut = true
            for cut2 in new_cuts[cut_type]
                if cut1[2] == cut2[2]
                    unique_cut = false
                    break
                end
            end
            if unique_cut
                push!(new_cuts[cut_type], cut1)
            end
        end
    end

    # Base.Main.@infiltrate
    return new_cuts
end

function _generate_cuts_from_dual(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    # Shall we generate cuts?
    gen_feas_cuts = has_attribute_type(model, Benders.AbstractFeasibilityCutType)
    gen_opt_cuts = has_attribute_type(model, Benders.AbstractOptimalityCutType)
    gen_feas_cuts || gen_opt_cuts || return new_cuts

    cur_sol_main = model.info[:results][:main][:sol]::JuMP.Containers.DenseAxisArray

    for i in 1:(length(model.models) - 1)
        m_sub = Benders.sub(model; index=i)
        get(m_sub.ext, :dualization_is_dualized, false) || continue

        # TODO: better checks!
        if gen_feas_cuts && JuMP.primal_status(m_sub) == MOI.INFEASIBILITY_CERTIFICATE
            # Primal unbounded => construct a feasibility cut from the ray.
            exp_cut = JuMP.AffExpr(JuMP.value(m_sub[:obj_base]))
            for elem in m_sub.ext[:dualization_obj_param]
                vi = m_sub.ext[:dualization_var_to_vi][elem[1]]
                λ = JuMP.value(elem[2])
                JuMP.add_to_expression!(exp_cut, λ, Benders.main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end
            push!(new_cuts[:feasibility], (i, exp_cut))
        elseif gen_opt_cuts && JuMP.termination_status(m_sub) == MOI.OPTIMAL && JuMP.primal_status(m_sub) == MOI.FEASIBLE_POINT
            exp_cut = JuMP.AffExpr(
                JuMP.dual_status(m_sub) == MOI.FEASIBLE_POINT ? JuMP.dual_objective_value(m_sub) : JuMP.objective_value(m_sub)
            )
            for elem in m_sub.ext[:dualization_obj_param]
                vi = m_sub.ext[:dualization_var_to_vi][elem[1]]
                λ = JuMP.value(elem[2])
                JuMP.add_to_expression!(exp_cut, λ, Benders.main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end
            push!(new_cuts[:optimality], (i, exp_cut))
        else
            @warn "Could not create the requested cut type" sub_model = i gen_feas_cuts gen_opt_cuts
        end
    end

    return new_cuts
end

function _generate_cuts_from_primal(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    # Shall we generate cuts?
    gen_feas_cuts = has_attribute_type(model, Benders.AbstractFeasibilityCutType)
    gen_opt_cuts = has_attribute_type(model, Benders.AbstractOptimalityCutType)
    gen_feas_cuts || gen_opt_cuts || return new_cuts

    vis_main = model.vis[1]
    cur_sol_main = model.info[:results][:main][:sol]::JuMP.Containers.DenseAxisArray

    for i in 1:(length(model.models) - 1)
        m_sub = Benders.sub(model; index=i)
        get(m_sub.ext, :dualization_is_dualized, false) && continue

        vis_sub = model.vis[1 + i]
        y_sub = m_sub[:y]

        cut_type = Benders.check_cut_type(m_sub)

        if cut_type == :optimality && gen_opt_cuts
            exp_cut = JuMP.AffExpr(0.0)
            JuMP.add_to_expression!(exp_cut, model.info[:results][:subs][i][:obj_lb])

            # TODO: check if this works as expected
            # if has_attribute_type(model, BD_MainObjectiveCon)
            #     JuMP.add_to_expression!(exp_cut, Benders.main(model)[:obj])
            # end

            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.reduced_cost(y_sub[vi])
                JuMP.add_to_expression!(exp_cut, λ, Benders.main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:optimality], (i, exp_cut))
        elseif cut_type == :feasibility && gen_feas_cuts
            exp_cut = JuMP.AffExpr(model.info[:results][:subs][i][:obj_dual])  # TODO: what's the best way to access this?
            
            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.dual(JuMP.FixRef(y_sub[vi]))
                JuMP.add_to_expression!(exp_cut, λ, Benders.main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:feasibility], (i, exp_cut))
        else
            @warn "Could not create the requested cut type" sub_model = i cut_type gen_feas_cuts gen_opt_cuts
        end
    end

    return new_cuts
end

function generate_cuts(model::Benders.DecomposedModel)
    new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    # Shall we generate cuts?
    gen_feas_cuts = has_attribute_type(model, Benders.AbstractFeasibilityCutType)
    gen_opt_cuts = has_attribute_type(model, Benders.AbstractOptimalityCutType)
    gen_feas_cuts || gen_opt_cuts || return new_cuts

    vis_main = model.vis[1]
    cur_sol_main = model.info[:results][:main][:sol]::JuMP.Containers.DenseAxisArray

    for i in 1:(length(model.models) - 1)
        vis_sub = model.vis[1 + i]
        m_sub = Benders.sub(model; index=i)
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

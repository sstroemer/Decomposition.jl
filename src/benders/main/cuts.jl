function modify(model::Benders.DecomposedModel, attribute::Benders.Main.AbstractCutType)
    # If it is a feasibility cut type, check if we can actually build those later on.
    if attribute isa Benders.Main.AbstractFeasibilityCutType
        if !has_attribute_type(model, Solver.ExtractDualRay) || !get_attribute(model, Solver.ExtractDualRay).activate
            @error "Adding feasibility cuts requires the solver to be able to extract dual rays"
            return nothing
        end
    end

    # For optimality cuts, we need to setup the θ variable.
    if attribute isa Benders.Main.AbstractOptimalityCutType && !haskey(main(model), :θ)
        # First time configuring the optimality cut type => setup the θ variable.
        JuMP.@variable(main(model), θ[i = 1:(length(model.models) - 1)])
        JuMP.set_lower_bound.(θ, -1e4)           # TODO
    end

    add_attribute!(model, attribute)
    return nothing
end

function generate_cuts(model::Benders.DecomposedModel)
    new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    # Shall we generate cuts?
    gen_feas_cuts = has_attribute_type(model, Benders.Main.AbstractFeasibilityCutType)
    gen_opt_cuts = has_attribute_type(model, Benders.Main.AbstractOptimalityCutType)
    gen_feas_cuts || gen_opt_cuts || return new_cuts

    vis_main = model.vis[1]
    cur_sol_main = model.info[:results][:main][:sol]::JuMP.Containers.DenseAxisArray

    for i in 1:(length(model.models) - 1)
        vis_sub = model.vis[1 + i]
        m_sub = sub(model; index=i)
        x_sub = m_sub[:x]

        cut_type = Benders.check_cut_type(m_sub)

        if cut_type == :optimality && gen_opt_cuts
            exp_cut = JuMP.AffExpr(0.0)
            JuMP.add_to_expression!(exp_cut, model.info[:results][:subs][i][:obj_lb])

            # TODO: check if this works as expected
            # if has_attribute_type(model, BD_MainObjectiveCon)
            #     JuMP.add_to_expression!(exp_cut, main(model)[:obj])
            # end

            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.reduced_cost(x_sub[vi])
                JuMP.add_to_expression!(exp_cut, λ, main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:optimality], (i, exp_cut))
        elseif cut_type == :feasibility && gen_feas_cuts
            exp_cut = JuMP.AffExpr(model.info[:results][:subs][i][:obj_dual])  # TODO: what's the best way to access this?
            
            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.dual(JuMP.FixRef(x_sub[vi]))
                JuMP.add_to_expression!(exp_cut, λ, main(model)[:x][vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:feasibility], (i, exp_cut))
        else
            @warn "Could not create the requested cut type" sub_model = i cut_type gen_feas_cuts gen_opt_cuts
        end
    end

    return new_cuts
end

function add_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    if !has_attribute_type(model, Benders.Main.AbstractCutType)
        @error "No cut type(s) specified"
        return 0, 0
    end

    nof_added_feas_cuts = 0
    nof_added_opt_cuts = 0

    if !isempty(new_cuts[:feasibility])
        # Use the "last", the currently active, cut type specfication.
        nof_added_feas_cuts = _add_feasibility_cuts(model, new_cuts, get_attributes(model, Benders.Main.AbstractFeasibilityCutType)[end])
    end

    if !isempty(new_cuts[:optimality])
        # Use the "last", the currently active, cut type specfication.
        nof_added_opt_cuts = _add_optimality_cuts(model, new_cuts, get_attributes(model, Benders.Main.AbstractOptimalityCutType)[end])
    end

    return nof_added_feas_cuts, nof_added_opt_cuts
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.FeasibilityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:feasibility]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(model.cuts[:feasibility], JuMP.@constraint(main(model), exp_agg_cut <= 0))
    
    return 1
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.FeasibilityCutTypeMulti)
    for (_, exp_cut) in new_cuts[:feasibility]
        push!(model.cuts[:feasibility], JuMP.@constraint(main(model), exp_cut <= 0))
    end

    return length(new_cuts[:feasibility])
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.FeasibilityCutTypeAggregated)
    @error "FeasibilityCutTypeAggregated not implemented"
    return 0
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.FeasibilityCutTypeAdaptive)
    @error "FeasibilityCutTypeAdaptive not implemented"
    return 0
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.OptimalityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:optimality]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(model.cuts[:optimality], JuMP.@constraint(main(model), sum(main(model)[:θ]) >= exp_agg_cut))
    
    return 1
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.OptimalityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:optimality]
        push!(model.cuts[:optimality], JuMP.@constraint(main(model), main(model)[:θ][i] >= exp_cut))
    end

    return length(new_cuts[:optimality])
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.OptimalityCutTypeAggregated)
    # TODO: it can happen that we add a lot of `\theta >= 0` cuts in OptimalityCutTypeMulti, due to the "trivial" sub-models.
    #       => catch linear dependent cuts and do not add them.
    # TODO: aggregating cuts is basically "averaging" them; we could: (1) weight the average, (2) cluster them before aggregating to generate `N >= 1` cuts
    @error "OptimalityCutTypeAggregated not implemented"
    return 0
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.OptimalityCutTypeAdaptive)
    @error "OptimalityCutTypeAdaptive not implemented"
    return 0
end

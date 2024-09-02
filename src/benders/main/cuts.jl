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

function _preprocess_cuts_remove_redundant!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    !has_attribute_type(model, Benders.CutPreprocessingRemoveRedundant) && return nothing

    # TODO: feasibility cuts can also be redundant within an iteration/ between sub-models: (1, -x[2905] + 31965.28), (4, -x[2905] + 36838.016), ...
    for type in [:feasibility, :optimality]
        valid = [true for _ in new_cuts[type]]

        for i in eachindex(new_cuts[type])
            cut = new_cuts[type][i]
            
            iter_oldcuts = (c for c in model.cuts[type] if c.sub_model_index == cut[1])
            isempty(iter_oldcuts) && continue

            last_opt_cut = last(c for c in model.cuts[type] if c.sub_model_index == cut[1])
            cut[2] != last_opt_cut.cut_exp && continue

            # This cut is identical to an old one.
            valid[i] = false
        end

        new_cuts[type] = [new_cuts[type][i] for i in eachindex(new_cuts[type]) if valid[i]]
    end

    return nothing
end

function _preprocess_cuts_stabilize_numerical_range!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    !has_attribute_type(model, Benders.CutPreprocessingStabilizeNumericalRange) && return nothing
    attr = get_attribute(model, Benders.CutPreprocessingStabilizeNumericalRange)

    # Remove extremely low coefficients.
    # Example: (4, -0.0016145468734458735 x[121] - 765000.7505230922 x[363] ... + 2.5318715674811035e10)

    for cut in new_cuts[:optimality]
        constant = cut[2].constant
        constant_factor = abs(constant) ./ abs.(cut[2].terms.vals)

        for (var, coeff, fact) in zip(cut[2].terms.keys, cut[2].terms.vals, constant_factor)
            if fact > attr.const_factor_threshold
                max_delta = Inf
                
                if coeff > 0 && JuMP.has_lower_bound(var)
                    max_delta = JuMP.lower_bound(var) * coeff
                elseif coeff < 0 && JuMP.has_upper_bound(var)
                    max_delta = JuMP.upper_bound(var) * coeff
                end
                
                if abs(max_delta / constant) < attr.const_factor_elimination_max_rel_delta
                    delete!(cut[2].terms, var)
                    cut[2].constant += max_delta
                end
            end
        end
    end

    return nothing
end

function preprocess_cuts!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    # TODO: keep stats on number of preprocessing steps, etc.
    _preprocess_cuts_remove_redundant!(model, new_cuts)
    _preprocess_cuts_stabilize_numerical_range!(model, new_cuts)

    return nothing
end

function postprocess_cuts!(model::Benders.DecomposedModel)
    return nothing
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
    
    push!(
        model.cuts[:feasibility],
        Benders.GeneralFeasibilityCut(
            iteration = current_iteration(model),
            cut_exp = exp_agg_cut,
            cut_con = JuMP.@constraint(main(model), exp_agg_cut <= 0)
        )
    )
    
    return 1
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.FeasibilityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:feasibility]
        push!(
            model.cuts[:feasibility],
            Benders.GeneralFeasibilityCut(
                iteration = current_iteration(model),
                sub_model_index = i,
                cut_exp = exp_cut,
                cut_con = JuMP.@constraint(main(model), exp_cut <= 0)
            )
        )
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
    
    push!(
        model.cuts[:optimality],
        Benders.GeneralOptimalityCut(
            iteration = current_iteration(model),
            cut_exp = exp_agg_cut,
            cut_con = JuMP.@constraint(main(model), sum(main(model)[:θ]) >= exp_agg_cut)
        )
    )
    
    return 1
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.Main.OptimalityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:optimality]
        push!(
            model.cuts[:optimality],
            Benders.GeneralOptimalityCut(
                iteration = current_iteration(model),
                sub_model_index = i,
                cut_exp = exp_cut,
                cut_con = JuMP.@constraint(main(model), main(model)[:θ][i] >= exp_cut)
            )
        )
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

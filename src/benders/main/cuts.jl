abstract type AbstractFeasibilityCutType <: AbstractCutType end
abstract type AbstractOptimalityCutType <: AbstractCutType end

@kwdef struct FeasibilityCutTypeSingle <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeMulti <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeAggregated <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeAdaptive <: AbstractFeasibilityCutType; end

@kwdef struct OptimalityCutTypeSingle <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeMulti <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeAggregated <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeAdaptive <: AbstractOptimalityCutType; end

function modify(model::DecomposedModel, attribute::AbstractCutType)
    setup_θ = !has_attribute_type(model, typeof(attribute))

    add_attribute!(model, attribute)

    if setup_θ
        # First time configuring the cut type => setup the θ variable.
        JuMP.@variable(main(model), θ[i = 1:(length(model.models) - 1)])
        JuMP.set_lower_bound.(θ, -1e4)           # TODO
    end

    return nothing
end

function generate_cuts(model::DecomposedModel)
    cur_sol_main = model.solutions[:current][:main]::JuMP.Containers.DenseAxisArray

    new_cuts = Dict{Symbol, Vector{Any}}(
        :feasibility => [],
        :optimality => []
    )

    x_main = main(model)[:x]
    vis_main = model.vis[1]

    for i in 1:(length(model.models) - 1)
        vis_sub = model.vis[1 + i]
        m_sub = sub(model; index=i)
        x_sub = m_sub[:x]

        cut_type = check_cut_type(m_sub)

        if cut_type == :optimality
            exp_cut = JuMP.AffExpr(0.0)
            JuMP.add_to_expression!(exp_cut, model.info[:results][:subs][i][:obj])  # TODO: this should be the "lower bound of the sub-model"

            # TODO: check if this works as expected
            # if has_attribute_type(model, BD_MainObjectiveCon)
            #     JuMP.add_to_expression!(exp_cut, main(model)[:obj])
            # end

            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.reduced_cost(x_sub[vi])
                JuMP.add_to_expression!(exp_cut, λ, main(model)[vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:optimality], (i, exp_cut))
        elseif cut_type == :feasibility
            exp_cut = JuMP.AffExpr(model.info[:results][:subs][i][:obj_dual])  # TODO: what's the best way to access this?
            
            for vi in vis_main
                (vi in vis_sub) || continue

                λ = JuMP.dual(JuMP.FixRef(x_sub[vi]))
                JuMP.add_to_expression!(exp_cut, λ, main(model)[vi])
                JuMP.add_to_expression!(exp_cut, λ, -cur_sol_main[vi])
            end

            push!(new_cuts[:feasibility], (i, exp_cut))
        else
            @error "Cut type for sub-model $(i) could not be determined"
        end
    end

    return new_cuts
end

function add_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    if !has_attribute_type(model, AbstractCutType)
        @error "No cut type(s) specified"
        return 0, 0
    end

    nof_added_feas_cuts = 0
    nof_added_opt_cuts = 0

    if !isempty(new_cuts[:feasibility]) && has_attribute_type(model, AbstractFeasibilityCutType)
        # Use the "last", the currently active, cut type specfication.
        nof_added_feas_cuts = _add_feasibility_cuts(model, new_cuts, get_attributes(model, AbstractFeasibilityCutType)[end])
    end

    if !isempty(new_cuts[:optimality]) && has_attribute_type(model, AbstractOptimalityCutType)
        # Use the "last", the currently active, cut type specfication.
        nof_added_opt_cuts = _add_optimality_cuts(model, new_cuts, get_attributes(model, AbstractOptimalityCutType)[end])
    end

    return nof_added_feas_cuts, nof_added_opt_cuts
end

function _add_feasibility_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::FeasibilityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:feasibility]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(model.cuts[:feasibility], JuMP.@constraint(main(model), exp_agg_cut <= 0))
    
    return 1
end

function _add_feasibility_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::FeasibilityCutTypeMulti)
    for (_, exp_cut) in new_cuts[:feasibility]
        push!(model.cuts[:feasibility], JuMP.@constraint(main(model), exp_cut <= 0))
    end

    return length(new_cuts[:feasibility])
end

function _add_feasibility_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::FeasibilityCutTypeAggregated)
    @error "FeasibilityCutTypeAggregated not implemented"
    return 0
end

function _add_feasibility_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::FeasibilityCutTypeAdaptive)
    @error "FeasibilityCutTypeAdaptive not implemented"
    return 0
end

function _add_optimality_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::OptimalityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:optimality]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(model.cuts[:optimality], JuMP.@constraint(main(model), sum(main(model)[:θ]) >= exp_agg_cut))
    
    return 1
end

function _add_optimality_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::OptimalityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:optimality]
        push!(model.cuts[:optimality], JuMP.@constraint(main(model), main(model)[:θ][i] >= exp_cut))
    end

    return length(new_cuts[:optimality])
end

function _add_optimality_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::OptimalityCutTypeAggregated)
    # TODO: it can happen that we add a lot of `\theta >= 0` cuts in OptimalityCutTypeMulti, due to the "trivial" sub-models.
    #       => catch linear dependent cuts and do not add them.
    # TODO: aggregating cuts is basically "averaging" them; we could: (1) weight the average, (2) cluster them before aggregating to generate `N >= 1` cuts
    @error "OptimalityCutTypeAggregated not implemented"
    return 0
end

function _add_optimality_cuts(model::DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::OptimalityCutTypeAdaptive)
    @error "OptimalityCutTypeAdaptive not implemented"
    return 0
end

function check_cut_type(jump_model::JuMP.Model; verbose::Bool = true)
    if JuMP.is_solved_and_feasible(jump_model)
        return :optimality
    end

    if JuMP.dual_status(jump_model) != JuMP.MOI.INFEASIBILITY_CERTIFICATE
        verbose && (@error "Turn off presolve, or any setting blocking extraction of dual rays")
        return :error
    end

    return :feasibility
end

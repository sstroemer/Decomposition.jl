function apply!(model::Benders.DecomposedModel, attribute::Benders.AbstractCutType)
    # If it is a feasibility cut type, check if we can actually build those later on.
    if attribute isa Benders.AbstractFeasibilityCutType
        if !has_attribute_type(model, Solver.ExtractDualRay) || !get_attribute(model, Solver.ExtractDualRay).activate
            @error "Adding feasibility cuts requires the solver to be able to extract dual rays"
            return false
        end
    end

    jm = Benders.main(model)

    # For optimality cuts, we need to setup the θ variable.
    if attribute isa Benders.AbstractOptimalityCutType && !haskey(jm, :θ)
        # First time configuring the optimality cut type => setup the θ variable.
        JuMP.@variable(jm, θ[i = 1:(length(model.models) - 1)])
        JuMP.set_lower_bound.(θ, -1e4)           # TODO
    end

    return true
end

# Cut processing is done as "flag" and is never "applied" to the model.
apply!(model::Benders.DecomposedModel, attribute::Benders.AbstractCutProcessing) = true

# TODO: transform the following into execution queries

include("preprocessing.jl")
include("generation.jl")
include("postprocessing.jl")

function add_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    if !has_attribute_type(model, Benders.AbstractCutType)
        @error "No cut type(s) specified"
        return 0, 0
    end

    nof_added_feas_cuts = 0
    nof_added_opt_cuts = 0

    if !isempty(new_cuts[:feasibility])
        # Use the "last", the currently active, cut type specfication.
        nof_added_feas_cuts = _add_feasibility_cuts(model, new_cuts, get_attributes(model, Benders.AbstractFeasibilityCutType)[end])
    end

    if !isempty(new_cuts[:optimality])
        # Use the "last", the currently active, cut type specfication.
        nof_added_opt_cuts = _add_optimality_cuts(model, new_cuts, get_attributes(model, Benders.AbstractOptimalityCutType)[end])
    end

    return nof_added_feas_cuts, nof_added_opt_cuts
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.FeasibilityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:feasibility]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(
        model.cuts[:feasibility],
        Benders.GeneralFeasibilityCut(
            iteration = current_iteration(model),
            cut_exp = exp_agg_cut,
            cut_con = JuMP.@constraint(Benders.main(model), exp_agg_cut <= 0)
        )
    )
    
    return 1
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.FeasibilityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:feasibility]
        push!(
            model.cuts[:feasibility],
            Benders.GeneralFeasibilityCut(
                iteration = current_iteration(model),
                sub_model_index = i,
                cut_exp = exp_cut,
                cut_con = JuMP.@constraint(Benders.main(model), exp_cut <= 0)
            )
        )
    end

    return length(new_cuts[:feasibility])
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.FeasibilityCutTypeAggregated)
    @error "FeasibilityCutTypeAggregated not implemented"
    return 0
end

function _add_feasibility_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.FeasibilityCutTypeAdaptive)
    @error "FeasibilityCutTypeAdaptive not implemented"
    return 0
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.OptimalityCutTypeSingle)
    exp_agg_cut = JuMP.AffExpr(0.0)
    
    for (_, exp_cut) in new_cuts[:optimality]
        JuMP.add_to_expression!(exp_agg_cut, exp_cut)
    end
    
    push!(
        model.cuts[:optimality],
        Benders.GeneralOptimalityCut(
            iteration = current_iteration(model),
            cut_exp = exp_agg_cut,
            cut_con = JuMP.@constraint(Benders.main(model), sum(Benders.main(model)[:θ]) >= exp_agg_cut)
        )
    )
    
    return 1
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.OptimalityCutTypeMulti)
    for (i, exp_cut) in new_cuts[:optimality]
        push!(
            model.cuts[:optimality],
            Benders.GeneralOptimalityCut(
                iteration = current_iteration(model),
                sub_model_index = i,
                cut_exp = exp_cut,
                cut_con = JuMP.@constraint(Benders.main(model), Benders.main(model)[:θ][i] >= exp_cut)
            )
        )
    end

    return length(new_cuts[:optimality])
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.OptimalityCutTypeAggregated)
    # TODO: it can happen that we add a lot of `\theta >= 0` cuts in OptimalityCutTypeMulti, due to the "trivial" sub-models.
    #       => catch linear dependent cuts and do not add them.
    # TODO: aggregating cuts is basically "averaging" them; we could: (1) weight the average, (2) cluster them before aggregating to generate `N >= 1` cuts
    @error "OptimalityCutTypeAggregated not implemented"
    return 0
end

function _add_optimality_cuts(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}}, ::Benders.OptimalityCutTypeAdaptive)
    @error "OptimalityCutTypeAdaptive not implemented"
    return 0
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.AbstractCutType)
    # If it is a feasibility cut type, check if we can actually build those later on.
    if attribute isa Benders.AbstractFeasibilityCutType
        if !has_attribute_type(model, Solver.ExtractDualRay) || !get_attribute(model, Solver.ExtractDualRay).activate
            @error "Adding feasibility cuts requires the solver to be able to extract dual rays"
            return false
        end
    end

    # For optimality cuts, we need to setup the θ variable.
    if attribute isa Benders.AbstractOptimalityCutType && !haskey(Benders.main(model), :θ)
        # First time configuring the optimality cut type => setup the θ variable.
        JuMP.@variable(Benders.main(model), θ[i = 1:(length(model.models) - 1)])
        JuMP.set_lower_bound.(θ, -1e4)           # TODO
    end

    return true
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.CutTypeMISFSZ)
    return true # TODO
    # if has_attribute_type(model, Benders.AbstractFeasibilityCutType)
    #     @error "Active cut type `MISFSZ` cannot be used with general feasibility cuts"
    #     return false
    # end

    # if has_attribute_type(model, Benders.AbstractOptimalityCutType)
    #     @error "Active cut type `MISFSZ` cannot be used with general optimality cuts"
    #     return false
    # end

    # for jm in Benders.subs(model)
    #     if !get(jm.ext, :dualization_is_dualized, false)
    #         @error "Active cut type `MISFSZ` requires the sub-models to be dualized"
    #         return false
    #     end

    #     # TODO: Make that a parameter
    #     ω_0 = 1.0

    #     # Add the π_0 variable to the sub-model.
    #     JuMP.@variable(jm, π_0 >= 0)

    #     # Multiply `d` (RHS) with `π_0`.
    #     # TODO: is it safe to leave out the `x >= ...` and `x <= ...` constraints?
    #     cons = JuMP.all_constraints(jm; include_variable_in_set_constraints = false)
    #     for con in cons
    #         co = JuMP.constraint_object(con)

    #         rhs = (co.set isa MOI.EqualTo) ? co.set.value : ((co.set isa MOI.LessThan) ? co.set.upper : co.set.lower)
    #         # JuMP.set_normalized_coefficient(con, π_0, -rhs)
    #         # JuMP.set_normalized_rhs(con, 0)           
    #     end

    #     # Add the "CGLP normalization condition" (the "reduced" version as suggested).
    #     JuMP.@constraint(
    #         jm,
    #         cglp_normalization,
    #         sum(elem[2] for elem in jm.ext[:dualization_obj_param]) + ω_0 * π_0 == 1
    #     )

    #     # Construct the additional objective term.
    #     jm[:obj_misfsz] = JuMP.AffExpr(0.0)
    # end

    # # Finally, we need the standard θ variables in the main-model (if not already constructed).
    # if !haskey(Benders.main(model), :θ)
    #     JuMP.@variable(Benders.main(model), θ[i = eachindex(Benders.subs(model))])
    #     JuMP.set_lower_bound.(θ, -1e4)           # TODO
    # end

    # return true
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

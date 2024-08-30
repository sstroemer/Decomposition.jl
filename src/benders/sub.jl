struct BD_SubObjectivePure <: DecompositionAttribute
    index::Int64

    BD_SubObjectivePure() = new(-1)
end

struct BD_SubObjectiveFull <: DecompositionAttribute
    index::Int64

    BD_SubObjectiveFull() = new(-1)
end

struct BD_SubEnsureFeasibilityFull <: DecompositionAttribute
    index::Int64
    penalty::Float64
    penalty_map::Dict{ConstraintRef, AffExpr}

    BD_SubEnsureFeasibility(index::Int64, penalty::Float64) = new(index, penalty, Dict{ConstraintRef, AffExpr}())
    BD_SubEnsureFeasibility(index::Int64) = new(index, 1e8, Dict{ConstraintRef, AffExpr}())
end

struct BD_SubEnsureFeasibilityLinked <: DecompositionAttribute
    index::Int64
    penalty::Float64

    BD_SubEnsureFeasibilityLinked(index::Int64, penalty::Float64) = new(index, penalty)
    BD_SubEnsureFeasibilityLinked() = new(-1, 1e8)
end

struct BD_SubEnsureFeasibilityRegex <: DecompositionAttribute
    index::Int64
    penalty::Float64
    regex::Regex
end

struct BD_SubFeasiblityCutsAlways <: DecompositionAttribute end
struct BD_SubFeasiblityCutsOnDemand <: DecompositionAttribute
    init_on::Bool
    state::Vector{Bool}

    BD_SubFeasiblityCutsOnDemand(init_on::Bool) = new(init_on, [])
end

function bd_jump_configure_dualrays(model::DecomposedModel, sub_model_index::Int, activate::Bool)
    @debug "Configure dual ray extraction" sub_model_index activate

    jump_model = bd_sub(model; index=sub_model_index)
    
    solver = solver_name(jump_model)
    if solver == "Gurobi"
        set_attribute(jump_model, "InfUnbdInfo", activate ? 1 : 0)
        set_attribute(jump_model, "DualReductions", activate ? 0 : 1)
    elseif solver == "HiGHS"
        set_attribute(jump_model, "presolve", activate ? "off" : "on")
        # set_attribute(jump_model, "solver", activate ? "simplex" : "choose")
        # TODO: does that also need "sovler=simplex" and "simplex_strategy=dual"?
        # See:  https://ampl.com/colab/notebooks/ampl-development-tutorial-56-parallelizing-subproblem-solves-in-benders-decomposition.html#utility-and-worker-functions
    else
        @error "Controlling dual ray extraction via `JuMP` is currently not supported for solver `$(solver)`"
        return nothing
    end

    if bd_has_attribute_type(model, BD_SubFeasiblityCutsOnDemand)
        attr = bd_get_attribute(model, BD_SubFeasiblityCutsOnDemand)
        attr.state[sub_model_index] = activate
    end

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubFeasiblityCutsAlways)
    # We need to push first, because `bd_jump_configure_dualrays` will look for this attribute.
    push!(model.attributes, attribute)

    for i in eachindex(bd_subs(model))
        bd_jump_configure_dualrays(model, i, true)
    end

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubFeasiblityCutsOnDemand)
    # We need to push first, because `bd_jump_configure_dualrays` will look for this attribute.
    push!(model.attributes, attribute)

    for i in eachindex(bd_subs(model))
        push!(attribute.state, false)
        bd_jump_configure_dualrays(model, i, attribute.init_on)
    end

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectivePure)
    vis_main = model.idx_model_vars[1]
    for i in 1:(length(model.models) - 1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = bd_sub(model; index=i)
        vis_sub = model.idx_model_vars[1 + i]
    
        obj = AffExpr(0.0)
        for vi in vis_sub
            (vi in vis_main) && continue
            add_to_expression!(obj, model.lpmd.c[vi], m_sub[:x][vi])
        end
        m_sub[:obj] = obj
        @objective(m_sub, Min, obj)
    end

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectiveFull)
    @error "This does not work yet (multiple subs need to correctly SHARE main-variable costs)!"

    m = bd_sub(model; index=attribute.index)
    idx_v = model.idx_model_vars[1 + attribute.index]   # TODO: account for "share of main variables, with multiple sub-models"

    obj = AffExpr(0.0)
    for i in idx_v
        # (i in model.idx_model_vars[1]) && continue    # TODO: divide here!
        add_to_expression!(obj, model.lpmd.c[i], m[:x][i])
    end
    @objective(m, Min, obj)

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubEnsureFeasibilityFull)
    @warn "Untested ..."

    for i in 1:(length(model.models) - 1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = bd_sub(model; index=i)
        relax_with_penalty!(m_sub; default = attribute.penalty)
    end

    # TODO: consider everything below ...

    # Get all variables that are part of the main-model, and find the constraints related to them.
    # idx_v = model.idx_model_vars[1]
    # idx_c_torelax = findall(sum(model.lpmd.A[:, idx_v] .!= 0; dims=2)[:, 1] .!= 0)

    # TODO: only relax the necessary constraints ...
    # relax_with_penalty!(bd_sub(model), Dict(linking_constraints .=> penalty))
    # merge!(attribute.penalty_map, relax_with_penalty!(m; default = attribute.penalty))
    
    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubEnsureFeasibilityLinked)
    vis_main = model.idx_model_vars[1]

    for i in 1:(length(model.models) - 1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = bd_sub(model; index=i)
        vis_main_in_sub = [vi for vi in vis_main if vi in m_sub[:x].axes[1]]

        @variable(m_sub, z_pos[i = vis_main_in_sub], lower_bound = 0)
        @variable(m_sub, z_neg[i = vis_main_in_sub], lower_bound = 0)

        all_con = all_constraints(m_sub; include_variable_in_set_constraints = false)

        exp_penalty = AffExpr(0.0)   
        for vi in vis_main_in_sub
            coeffs = normalized_coefficient.(all_con, m_sub[:x][vi])
            set_normalized_coefficient.(all_con, z_pos[vi], coeffs)
            set_normalized_coefficient.(all_con, z_neg[vi], -coeffs)

            add_to_expression!(exp_penalty, z_pos[vi], attribute.penalty)
            add_to_expression!(exp_penalty, z_neg[vi], attribute.penalty)
        end

        @objective(m_sub, Min, objective_function(m_sub) + exp_penalty)
    end

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubEnsureFeasibilityRegex)
    vis_main = model.idx_model_vars[1]
    vn_main = name.(model.lpmd.variables[vis_main])
    vis_main = [vi for (vi, vn) in zip(vis_main, vn_main) if !isnothing(match(attribute.regex, vn))]

    for i in 1:(length(model.models) - 1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = bd_sub(model; index=i)
        vis_main_in_sub = [vi for vi in vis_main if vi in m_sub[:x].axes[1]]

        @variable(m_sub, z_pos[i = vis_main_in_sub], lower_bound = 0)
        @variable(m_sub, z_neg[i = vis_main_in_sub], lower_bound = 0)

        all_con = all_constraints(m_sub; include_variable_in_set_constraints = false)

        exp_penalty = AffExpr(0.0)   
        for vi in vis_main_in_sub
            coeffs = normalized_coefficient.(all_con, m_sub[:x][vi])
            set_normalized_coefficient.(all_con, z_pos[vi], coeffs)
            set_normalized_coefficient.(all_con, z_neg[vi], -coeffs)

            add_to_expression!(exp_penalty, z_pos[vi], attribute.penalty)
            add_to_expression!(exp_penalty, z_neg[vi], attribute.penalty)
        end

        @objective(m_sub, Min, objective_function(m_sub) + exp_penalty)
    end

    push!(model.attributes, attribute)
    return nothing
end


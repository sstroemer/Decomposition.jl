struct BD_SubObjectivePure <: DecompositionAttribute; index::Int64; end
struct BD_SubObjectiveFull <: DecompositionAttribute; index::Int64; end

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
    BD_SubEnsureFeasibilityLinked(index::Int64) = new(index, 1e8)
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectivePure)
    m = bd_sub(model; index=attribute.index)
    idx_v = model.idx_model_vars[1 + attribute.index]

    obj = AffExpr(0.0)
    for i in idx_v
        (i in model.idx_model_vars[1]) && continue
        add_to_expression!(obj, model.lpmd.c[i], m[:x][i])
    end
    @objective(m, Min, obj)

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectiveFull)
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
    m = bd_sub(model; index=attribute.index)

    # Get all variables that are part of the main-model, and find the constraints related to them.
    # idx_v = model.idx_model_vars[1]
    # idx_c_torelax = findall(sum(model.lpmd.A[:, idx_v] .!= 0; dims=2)[:, 1] .!= 0)

    # TODO: only relax the necessary constraints ...
    # relax_with_penalty!(bd_sub(model), Dict(linking_constraints .=> penalty))
    merge!(attribute.penalty_map, relax_with_penalty!(m; default = attribute.penalty))
    
    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubEnsureFeasibilityLinked)
    m = bd_sub(model; index=attribute.index)

    # Get all variables that are part of the main-model, and find the constraints related to them.
    idx_v = model.idx_model_vars[1]

    @variable(m, z_pos[i = idx_v], lower_bound = 0)
    @variable(m, z_neg[i = idx_v], lower_bound = 0)

    all_con = all_constraints(m; include_variable_in_set_constraints = false)

    exp_penalty = AffExpr(0.0)   
    for i in idx_v
        coeffs = normalized_coefficient.(all_con, m[:x][i])
        set_normalized_coefficient.(all_con, z_pos[i], coeffs)
        set_normalized_coefficient.(all_con, z_neg[i], -coeffs)

        add_to_expression!(exp_penalty, z_pos[i], attribute.penalty)
        add_to_expression!(exp_penalty, z_neg[i], attribute.penalty)
    end

    @objective(m, Min, objective_function(m) + exp_penalty)

    push!(model.attributes, attribute)
    return nothing
end

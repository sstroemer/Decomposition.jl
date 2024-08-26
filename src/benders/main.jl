struct BD_MainObjectiveConstOnly <: DecompositionAttribute end
struct BD_MainObjectiveObj <: DecompositionAttribute end
struct BD_MainObjectiveCon <: DecompositionAttribute end

struct MainVirtualBounds <: DecompositionAttribute
    lower::Float64
    upper::Float64
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

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveConstOnly)
    m = bd_main(model)

    if !haskey(bd_main(model), :θ)
        @variable(bd_main(model), θ)
        set_lower_bound(θ, 0)  # TODO
    end
    
    @objective(m, Min, model.lpmd.c_offset + θ)

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveObj)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]

    if !haskey(bd_main(model), :θ)
        @variable(bd_main(model), θ)
        set_lower_bound(θ, 0)  # TODO
    end
    
    @objective(m, Min, model.lpmd.c[idx_v]' * m[:x].data + model.lpmd.c_offset + θ)

    push!(model.attributes, attribute)
    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_MainObjectiveCon)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]

    if !haskey(bd_main(model), :θ)
        @variable(bd_main(model), θ)
        set_lower_bound(θ, 0)  # TODO
    end
    
    @expression(m, obj, model.lpmd.c[idx_v]' * m[:x].data)
    @objective(m, Min, model.lpmd.c_offset + θ)

    push!(model.attributes, attribute)
    return nothing
end

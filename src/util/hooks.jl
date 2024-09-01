function add_hook_main(model::Benders.DecomposedModel, when::Symbol, attribute::AbstractDecompositionAttribute)
    return _add_hook(model, 1, when, attribute)
end

function add_hook_sub(model::Benders.DecomposedModel, index::Int, when::Symbol, attribute::AbstractDecompositionAttribute)
    return _add_hook(model, index + 1, when, attribute)
end

function _add_hook(model::Benders.DecomposedModel, index::Int, when::Symbol, attribute::AbstractDecompositionAttribute)
    if !haskey(model.models[index].ext, :hooks)
        model.models[index].ext[:hooks] = Dict()
    end

    if !haskey(model.models[index].ext[:hooks], when)
        model.models[index].ext[:hooks][when] = AbstractDecompositionAttribute[]
    end
    
    push!(model.models[index].ext[:hooks][when], attribute)

    return nothing
end

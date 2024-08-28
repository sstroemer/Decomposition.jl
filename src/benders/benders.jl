bd_main(model::DecomposedModel) = model.models[1]
bd_sub(model::DecomposedModel; index::Int = 1) = model.models[1 + index]
bd_subs(model::DecomposedModel) = model.models[2:end]

bd_has_attribute_type(model::DecomposedModel, type::Type{T}) where T <: DecompositionAttribute = any(attr isa type for attr in model.attributes)

function bd_get_attributes(model::DecomposedModel, type::Type{T}) where T <: DecompositionAttribute
    attributes = DecompositionAttribute[]
    for attr in model.attributes
        attr isa type && push!(attributes, attr)
    end

    if isempty(attributes)
        @error "Model does not contain an attribute of type `$(T)`, check first with `bd_has_attribute_type(...)`"
    end

    return attributes
end

function bd_get_attribute(model::DecomposedModel, type::Type{T}) where T <: DecompositionAttribute
    attributes = bd_get_attributes(model, type)

    if length(attributes) > 1
        @error "Attributes of type `$(T)` not unique, use `bd_get_attributes(...)` instead"
        return nothing
    end

    return attributes[1]
end

function bd_modify(::DecomposedModel, ::DecompositionAttribute)
    @error "Not implemented"
end

include("main.jl")
include("sub.jl")
include("queries.jl")

function bd_update_fixed_variables(model::DecomposedModel, current_solution::JuMP.Containers.DenseAxisArray)
    for i in 1:(length(model.models) - 1)
        m_sub = bd_sub(model; index=i)
        for j in current_solution.axes[1]
            (j in model.idx_model_vars[1 + i]) || continue
            fix(m_sub[:x][j], current_solution[j]; force=true)
        end
    end

    return nothing
end


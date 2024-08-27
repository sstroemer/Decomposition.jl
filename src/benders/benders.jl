bd_main(model::DecomposedModel) = model.models[1]
bd_sub(model::DecomposedModel; index::Int = 1) = model.models[1 + index]
bd_subs(model::DecomposedModel) = model.models[2:end]

bd_has_attribute_type(model::DecomposedModel, type::Type{T}) where T <: DecompositionAttribute = any(attr isa type for attr in model.attributes)

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


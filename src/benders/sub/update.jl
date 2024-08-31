@timeit function update_fixed_variables(model::Benders.DecomposedModel, current_solution::JuMP.Containers.DenseAxisArray)
    for i in 1:(length(model.models) - 1)
        m_sub = sub(model; index=i)
        for j in current_solution.axes[1]
            (j in model.idx_model_vars[1 + i]) || continue
            JuMP.fix(m_sub[:x][j], current_solution[j]; force=true)
        end
    end

    return nothing
end

struct BD_SubFeasibility <: DecompositionAttribute; end
function bd_query(model::DecomposedModel, attribute::BD_SubFeasibility)
    # TODO: account for different previous ways to ensure feasibility (by checking the list of attributes/modifications)
    violation = sum(value.(bd_sub(model)[:z_pos])) + sum(value.(bd_sub(model)[:z_neg]))
    return isapprox(violation, 0.0; atol=1e-6)
end

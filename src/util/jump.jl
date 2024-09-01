function jump_objective_lb(jump_model::JuMP.Model)
    JuMP.is_solved_and_feasible(jump_model) || return missing

    obj_primal = try; JuMP.objective_value(jump_model); catch; +Inf; end
    obj_dual = try; JuMP.dual_objective_value(jump_model); catch; +Inf; end
    obj_fval = try; JuMP.value(JuMP.objective_function(jump_model)); catch; +Inf; end

    return min(obj_primal, obj_dual, obj_fval)
end

function jump_objective_ub(jump_model::JuMP.Model)
    JuMP.is_solved_and_feasible(jump_model) || return missing

    obj_primal = try; JuMP.objective_value(jump_model); catch; -Inf; end
    obj_dual = try; JuMP.dual_objective_value(jump_model); catch; -Inf; end
    obj_fval = try; JuMP.value(JuMP.objective_function(jump_model)); catch; -Inf; end

    return max(obj_primal, obj_dual, obj_fval)
end

function jump_safe_objective_value(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)

    return try; return JuMP.objective_value(jump_model); catch; missing; end
end

function jump_safe_dual_objective_value(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)
    
    return try; return JuMP.dual_objective_value(jump_model); catch; missing; end
end

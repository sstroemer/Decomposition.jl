function _fast_est_dual_obj_val(jump_model::JuMP.Model)
    get_b(set::JuMP.MOI.GreaterThan) = set.lower
    get_b(set::JuMP.MOI.LessThan) = set.upper
    get_b(set::JuMP.MOI.EqualTo) = set.value

    all_con = JuMP.all_constraints(jump_model; include_variable_in_set_constraints = true)

    π = JuMP.dual.(all_con)
    b = get_b.(JuMP.MOI.get.(jump_model, JuMP.MOI.ConstraintSet(), all_con))

    return JuMP.objective_sense(jump_model) == JuMP.MAX_SENSE ? -(π' * b) : (π' * b)
end

function jump_objective_all(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)

    obj_primal = try; JuMP.objective_value(jump_model); catch; missing; end
    obj_dual = (
        if JuMP.solver_name(jump_model) == "HiGHS"
            try; _fast_est_dual_obj_val(jump_model); catch; missing; end
        else
            try; JuMP.dual_objective_value(jump_model); catch; missing; end
        end
    )
    obj_fval = try; JuMP.value(JuMP.objective_function(jump_model)); catch; missing; end

    return (primal=obj_primal, dual=obj_dual, fval=obj_fval)
end

function jump_objective_lb(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)

    obj_primal = try; JuMP.objective_value(jump_model); catch; +Inf; end
    obj_dual = (
        if JuMP.solver_name(jump_model) == "HiGHS"
            try; _fast_est_dual_obj_val(jump_model); catch; +Inf; end
        else
            try; JuMP.dual_objective_value(jump_model); catch; +Inf; end
        end
    )
    obj_fval = try; JuMP.value(JuMP.objective_function(jump_model)); catch; +Inf; end

    return min(obj_primal, obj_dual, obj_fval)
end

function jump_objective_ub(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)

    obj_primal = try; JuMP.objective_value(jump_model); catch; -Inf; end
    obj_dual = (
        if JuMP.solver_name(jump_model) == "HiGHS"
            try; _fast_est_dual_obj_val(jump_model); catch; +Inf; end
        else
            try; JuMP.dual_objective_value(jump_model); catch; +Inf; end
        end
    )
    obj_fval = try; JuMP.value(JuMP.objective_function(jump_model)); catch; -Inf; end

    return max(obj_primal, obj_dual, obj_fval)
end

function jump_safe_objective_value(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)

    return try; return JuMP.objective_value(jump_model); catch; missing; end
end

function jump_safe_dual_objective_value(jump_model::JuMP.Model; require_feasibility::Bool=true)
    require_feasibility && (JuMP.is_solved_and_feasible(jump_model) || return missing)
    
    return (
        if JuMP.solver_name(jump_model) == "HiGHS"
            try; _fast_est_dual_obj_val(jump_model); catch; +Inf; end
        else
            try; JuMP.dual_objective_value(jump_model); catch; +Inf; end
        end
    )
end

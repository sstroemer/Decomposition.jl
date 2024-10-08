function _gurobi_modify_jump(jump_model::JuMP.Model, attribute::AlgorithmSimplex)
    if attribute.mode == :primal
        JuMP.set_attribute(jump_model, "Method", 0)
        return true
    elseif attribute.mode == :dual
        JuMP.set_attribute(jump_model, "Method", 1)
        return true
    end

    return false
end

function _gurobi_modify_jump(jump_model::JuMP.Model, attribute::AlgorithmIPM)
    JuMP.set_attribute(jump_model, "Method", 2)
    JuMP.set_attribute(jump_model, "Crossover", attribute.crossover ? -1 : 0)
    return true
end

function _gurobi_modify_jump(jump_model::JuMP.Model, attribute::ExtractDualRay)
    JuMP.set_attribute(jump_model, "InfUnbdInfo", attribute.activate ? 1 : 0)
    JuMP.set_attribute(jump_model, "DualReductions", attribute.activate ? 0 : 1)

    if attribute.activate && JuMP.get_attribute(jump_model, "Crossover") == 0
        @warn "`ExtractDualRay` requires crossover to be enabled, overwriting the previous setting"
        JuMP.set_attribute(jump_model, "Crossover", -1)
    end

    return true
end

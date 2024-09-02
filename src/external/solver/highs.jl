function _highs_modify_jump(model::JuMP.Model, attribute::AlgorithmSimplex)
    JuMP.set_attribute(model, "solver", "simplex")

    if attribute.mode == :primal
        JuMP.set_attribute(model, "simplex_strategy", 4)
        return true
    elseif attribute.mode == :dual
        # TODO: Check if we should go for `2` (dual-pami) instead
        JuMP.set_attribute(model, "simplex_strategy", 1)
        return true
    end

    return false
end

function _highs_modify_jump(model::JuMP.Model, attribute::AlgorithmIPM)
    JuMP.set_attribute(model, "solver", "ipm")
    JuMP.set_attribute(model, "run_crossover", attribute.crossover ? "on" : "off")
    return true
end

function _highs_modify_jump(jump_model::JuMP.Model, attribute::ExtractDualRay)
    JuMP.set_attribute(jump_model, "presolve", attribute.activate ? "off" : "on")

    if attribute.activate && JuMP.get_attribute(jump_model, "run_crossover") == "off"
        @warn "`ExtractDualRay` requires crossover to be enabled, overwriting the previous setting"
        JuMP.set_attribute(jump_model, "run_crossover", "on")
    end

    # TODO: does that also need "sovler=simplex" and "simplex_strategy=dual"?
    # See:  https://ampl.com/colab/notebooks/ampl-development-tutorial-56-parallelizing-subproblem-solves-in-benders-decomposition.html#utility-and-worker-functions
    # set_attribute(jump_model, "solver", activate ? "simplex" : "choose")

    return true
end

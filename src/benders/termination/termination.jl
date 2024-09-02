function modify(model::Benders.DecomposedModel, attribute::Benders.Termination.AbstractCriterion)
    add_attribute!(model, attribute)
    return nothing
end

function check_termination(model::Benders.DecomposedModel)
    function _inner(model::Benders.DecomposedModel)
        has_attribute_type(model, Benders.Termination.AbstractCriterion) || return false
        termination = get_attributes(model, Benders.Termination.AbstractCriterion)[end]

        if termination isa Benders.Termination.Stop
            iterations = current_iteration(model)
            seconds = TimerOutputs.tottime(model.timer) / 1e9
            opt_gap_abs = best_gap_abs(model)
            opt_gap_rel = best_gap_rel(model)
    
            if termination.iterations >= 0 && iterations >= termination.iterations
                @info "Termination criterion reached: [iterations]" iterations seconds opt_gap_abs opt_gap_rel
                return true
            end
    
            if termination.seconds >= 0 && seconds >= termination.seconds
                @info "Termination criterion reached: [seconds]" iterations seconds opt_gap_abs opt_gap_rel
                return true
            end
    
            if termination.opt_gap_abs >= 0 && opt_gap_abs <= termination.opt_gap_abs
                @info "Termination criterion reached: [opt_gap_abs]" iterations seconds opt_gap_abs opt_gap_rel
                return true
            end
    
            if termination.opt_gap_rel >= 0 && opt_gap_rel <= termination.opt_gap_rel
                @info "Termination criterion reached: [opt_gap_rel]" iterations seconds opt_gap_abs opt_gap_rel
                return true
            end
        end

        return false
    end

    terminate = _inner(model)

    if terminate
        @info "Model statistics" lb = best_lower_bound(model) ub = best_upper_bound(model) ncuts_feas = length(model.cuts[:feasibility]) ncuts_opt = length(model.cuts[:optimality])
    end

    return terminate
end

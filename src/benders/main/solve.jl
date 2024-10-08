function execute!(model::Benders.DecomposedModel, query::Benders.Query.SolveMain)
    jm = Benders.main(model)

    if has_attribute_type(model, Benders.Main.RegularizationLevelSet)
        ub = best_upper_bound(model)

        # TODO: make that an attribute
        # TODO / WRITING: this helps especially with HiGHS, and when solving the dualized subs, but also when solving the normal/primal one; it also slightly helps Gurobi for both cases
        if !isfinite(ub)
            ub = 1e12
        end

        if !isfinite(ub)
            # We have no valid upper bound, so nothing we can do right now.
            JuMP.optimize!(jm)

            model.info[:results][:main][:obj_lb] = jump_objective_lb(jm)
        else
            reg_attr = get_attribute(model, Benders.Main.RegularizationLevelSet)

            # Get the current lower bound.
            lb = begin
                reset_run_crossover = nothing
                if JuMP.solver_name(jm) == "HiGHS"
                    # TODO: remove this workaround
                    reset_run_crossover = JuMP.get_attribute(jm, "run_crossover")
                    JuMP.set_attribute(jm, "run_crossover", "on")
                end
                
                JuMP.@objective(jm, Min, jm[:obj_full])
                # TODO: only have the `reg_levelset_constraint` active if doing the actual level-set calculation
                #       instead of:
                JuMP.is_fixed(jm[:reg_levelset_L]) && JuMP.unfix(jm[:reg_levelset_L])
                JuMP.optimize!(jm)
                
                candidate = jump_objective_lb(jm; require_feasibility=false)
                if ismissing(candidate)
                    # TODO: fix this
                    # This is a workaround for HiGHS not finding a valid dual solution form the previous warm-start, but finding one by just re-running...
                    JuMP.optimize!(jm)
                    candidate = jump_objective_lb(jm; require_feasibility=false)
                    # TODO: afterwards remove the `require_feasibility` setting above too!
                end

                if JuMP.solver_name(jm) == "HiGHS" && !isfinite(candidate)
                    # TODO: fix this
                    # This happens with HiGHS if crossover is "imprecise", leading to a simplex cleanup, which then
                    # might indicate an unbounded problem (even if it is not).
                    # Try again, this time without crossover ... and just take the primal result.
                    JuMP.set_attribute(jm, "run_crossover", "off")
                    JuMP.optimize!(jm)
                    candidate = JuMP.objective_value(jm)
                end

                if JuMP.solver_name(jm) == "HiGHS" && !isnothing(reset_run_crossover)
                    # TODO: remove this workaround
                    JuMP.set_attribute(jm, "run_crossover", reset_run_crossover)
                end

                candidate
            end

            if !isfinite(lb)
                lb = -1e8  # TODO: make that an attribute
            end

            model.info[:results][:main][:obj_lb] = lb

            # Switch to level-set mode.
            JuMP.@objective(jm, Min, 0)  # TODO: is there a better way to switch to feasibility mode?

            alpha = reg_attr.alpha
            for _ in 1:reg_attr.safety_max_infeasible_resolve
                JuMP.fix(
                    jm[:reg_levelset_L],
                    alpha * ub + (1 - alpha) * lb;
                    force=true
                )

                # Re-optimize.
                JuMP.optimize!(jm)

                JuMP.is_solved_and_feasible(jm) && break
                alpha = min(1.0, alpha + reg_attr.infeasible_alpha_step)
            end

            if !JuMP.is_solved_and_feasible(jm)
                # We were not able to find a feasible solution.
                @error "Unable to find a feasible solution to main-model, possibly due to level-set"
            end
        end
    else
        JuMP.optimize!(jm)
    end

    return true
end

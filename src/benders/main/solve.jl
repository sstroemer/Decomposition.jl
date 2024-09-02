function solve_main(model::Benders.DecomposedModel)
    if has_attribute_type(model, Benders.Main.RegularizationLevelSet)
        ub = best_upper_bound(model)
        if !isfinite(ub)
            # We have no valid upper bound, so nothing we can do right now.
            JuMP.optimize!(main(model))

            model.info[:results][:main][:obj_lb] = jump_objective_lb(main(model))
        else
            reg_attr = get_attribute(model, Benders.Main.RegularizationLevelSet)

            # Get the current lower bound.
            lb = begin
                reset_run_crossover = nothing
                if JuMP.solver_name(main(model)) == "HiGHS"
                    # TODO: remove this workaround
                    reset_run_crossover = JuMP.get_attribute(main(model), "run_crossover")
                    JuMP.set_attribute(main(model), "run_crossover", "on")
                end
                
                JuMP.@objective(main(model), Min, main(model)[:obj_full])
                # TODO: only have the `reg_levelset_constraint` active if doing the actual level-set calculation
                #       instead of:
                JuMP.is_fixed(main(model)[:reg_levelset_L]) && JuMP.unfix(main(model)[:reg_levelset_L])
                JuMP.optimize!(main(model))
                
                candidate = jump_objective_lb(main(model); require_feasibility=false)
                if ismissing(candidate)
                    # TODO: fix this
                    # This is a workaround for HiGHS not finding a valid dual solution form the previous warm-start, but finding one by just re-running...
                    JuMP.optimize!(main(model))
                    candidate = jump_objective_lb(main(model); require_feasibility=false)
                    # TODO: afterwards remove the `require_feasibility` setting above too!
                end

                if JuMP.solver_name(main(model)) == "HiGHS" && !isnothing(reset_run_crossover)
                    # TODO: remove this workaround
                    JuMP.set_attribute(main(model), "run_crossover", reset_run_crossover)
                end

                candidate
            end
            model.info[:results][:main][:obj_lb] = lb
            
            # @show lb ub JuMP.dual_status(main(model)) JuMP.has_duals(main(model)) JuMP.dual_objective_value(main(model))

            # Switch to level-set mode.
            JuMP.@objective(main(model), Min, 0)  # TODO: is there a better way to switch to feasibility mode?

            alpha = reg_attr.alpha
            for _ in 1:reg_attr.safety_max_infeasible_resolve
                JuMP.fix(
                    main(model)[:reg_levelset_L],
                    alpha * ub + (1 - alpha) * lb;
                    force=true
                )

                # Re-optimize.
                JuMP.optimize!(main(model))

                JuMP.is_solved_and_feasible(main(model)) && break
                alpha = min(1.0, alpha + reg_attr.infeasible_alpha_step)
            end

            if !JuMP.is_solved_and_feasible(main(model))
                # We were not able to find a feasible solution.
                @error "Unable to find a feasible solution to main-model, possibly due to level-set"
            end
        end
    else
        JuMP.optimize!(main(model))
    end

    return nothing
end

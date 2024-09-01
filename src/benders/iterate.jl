function iterate!(model::Benders.DecomposedModel)   
    @timeit model.timer "main" begin
        @timeit model.timer "optimize" solve_main(model)
        @timeit model.timer "extract solution" extract_main(model)
    end

    # TODO: we could abort here if the gap is small enough (saving one final sub-model iteration) ?

    # Turn off dual ray extraction for sub-models that we expect to be feasible.
    # while !isempty(turn_on_again)
    #     i = pop!(turn_on_again)
    #     bd_jump_configure_dualrays(model, i, false)

    #     # TODO:
    #     # MOI.Utilities.reset_optimizer(bd_sub(model; index=i))   # only works for non-direct mode
    # end

    @timeit model.timer "sub" begin
        for i in 1:(length(model.models) - 1)
            @timeit model.timer "[$i]" begin
                m_sub = sub(model; index=i)

                # Fix the current main-model solution in the sub-model.
                @timeit model.timer "fix variables" begin
                    for j in model.info[:results][:main][:sol].axes[1]
                        (j in model.vis[1 + i]) || continue
                        JuMP.fix(m_sub[:x][j], model.info[:results][:main][:sol][j]; force=true)
                    end
                end

                # TODO: Implement "on demand" feasibility cuts again
                # See: https://github.com/sstroemer/Decomposition.jl/blob/df8f334fc1dc148e54b1638d83b0240d788627df/src/dev_draft.jl#L610
                # These can then be tagged as "re-optimize"

                # Solve the sub-model.
                @timeit model.timer "optimize" JuMP.optimize!(m_sub)

                @timeit model.timer "extract solution" begin
                    push!(
                        model.info[:results][:subs],
                        Dict(
                            :obj => jump_safe_objective_value(m_sub),
                            :obj_dual => jump_safe_dual_objective_value(m_sub; require_feasibility=false),
                            :obj_lb => jump_objective_lb(m_sub),
                            :obj_ub => jump_objective_ub(m_sub),
                        )
                    )
                    nothing
                end
            end
        end
    end

    added_cuts = @timeit model.timer "main" begin
        new_cuts = @timeit model.timer "cuts (generate)" generate_cuts(model)
        @timeit model.timer "cuts (add)" add_cuts(model, new_cuts)
    end
    
    # Pass added cuts and timings (in nanoseconds) to `next_iteration!`.
    next_iteration!(model, added_cuts)

    # if best_gap_rel(model) <= 1e-2
    #     @info "Reached relative gap of 1e-2" iterations = current_iteration(model) wall_time = round(total_wall_time(model) / 1e9; digits=2) cpu_time = round(total_cpu_time(model) / 1e9; digits=2)
    #     break
    # end

    terminate = check_termination(model)
    if terminate
        @info "Model statistics" lb = best_lower_bound(model) ub = best_upper_bound(model) ncuts_feas = length(model.cuts[:feasibility]) ncuts_opt = length(model.cuts[:optimality])
    end
    return terminate
end

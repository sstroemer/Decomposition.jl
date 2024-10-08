function execute!(model::Benders.DecomposedModel, query::Benders.Query.SolveSub)
    # TODO: Implement "on demand" feasibility cuts again
    # See: https://github.com/sstroemer/Decomposition.jl/blob/df8f334fc1dc148e54b1638d83b0240d788627df/src/dev_draft.jl#L610
    # These can then be tagged as "re-optimize"

    # Turn off dual ray extraction for sub-models that we expect to be feasible.
    # while !isempty(turn_on_again)
    #     i = pop!(turn_on_again)
    #     bd_jump_configure_dualrays(model, i, false)

    #     # TODO:
    #     # MOI.Utilities.reset_optimizer(bd_sub(model; index=i))   # only works for non-direct mode
    # end

    jm = Benders.sub(model; index=query.index)

    # Fix the current main-model solution in the sub-model.
    @timeit model.timer "fix variables" begin
        if get(jm.ext, :dualization_is_dualized, false)
            # Dualized model.
            Base.Main.@infiltrate

            # Construct all objective parts.
            jm[:obj_base] = jm.ext[:dualization_obj_base]

            jm[:obj_param] = JuMP.AffExpr(0.0)
            for elem in jm.ext[:dualization_obj_param]
                JuMP.add_to_expression!(
                    jm[:obj_param],
                    model.info[:results][:main][:sol][v2v_map[elem[1]]] * elem[3],
                    elem[2]
                )
            end

            jm[:obj] = jm[:obj_param] + jm[:obj_base]

            # Update objective.
            JuMP.@objective(jm, Min, jm[:obj])
        else
            # Standard (primal) model.
            for vi in jm[:y].axes[1]
                JuMP.fix(jm[:y][vi], model.info[:results][:main][:sol][vi]; force=true)
            end
        end
    end
    
    # Solve the sub's JuMP model.
    @timeit model.timer "optimize" JuMP.optimize!(jm)

    return nothing
end

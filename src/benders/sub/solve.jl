function solve_sub(model::Benders.DecomposedModel; index::Int64)
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

    jump_model = sub(model; index=index)

    # Fix the current main-model solution in the sub-model.
    @timeit model.timer "fix variables" begin
        for vi in model.info[:results][:main][:sol].axes[1]
            (vi in model.vis[1 + index]) || continue
            JuMP.fix(jump_model[:x][vi], model.info[:results][:main][:sol][vi]; force=true)
        end
    end
    
    # Solve the sub-model.
    @timeit model.timer "optimize" JuMP.optimize!(jump_model)

    return nothing
end

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

    # TODO: rework fixing to "==", which even helps with relaxing the main-sub-link

    jm = Benders.sub(model; index=query.index)

    # Fix the current main-model solution in the sub-model.
    @timeit model.timer "fix variables" begin
        for vi in model.info[:results][:main][:sol].axes[1]
            (vi in model.vis[1 + query.index]) || continue
            JuMP.fix(jm[:x][vi], model.info[:results][:main][:sol][vi]; force=true)
        end
    end
    
    # Solve the sub's JuMP model.
    @timeit model.timer "optimize" JuMP.optimize!(jm)

    return nothing
end

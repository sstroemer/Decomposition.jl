function execute!(model::Benders.DecomposedModel, query::Benders.Query.ExtractResultsSub)
    @timeit model.timer "extract solution" begin
        # TODO: require_feasibility=true
        isaf = JuMP.is_solved_and_feasible(Benders.sub(model; index=query.index))
        sol_obj = jump_objective_all(Benders.sub(model; index=query.index); require_feasibility=false)

        model.info[:results][:subs][query.index][:obj] = isaf ? sol_obj.primal : missing
        model.info[:results][:subs][query.index][:obj_dual] = sol_obj.dual
        model.info[:results][:subs][query.index][:obj_lb] = isaf ? minimum(skipmissing(sol_obj)) : missing
        model.info[:results][:subs][query.index][:obj_ub] = isaf ? maximum(skipmissing(sol_obj)) : missing

        nothing
    end

    return true
end

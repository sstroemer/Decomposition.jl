function extract_sub(model::Benders.DecomposedModel; index::Int64)
    @timeit model.timer "extract solution" begin
        # TODO: require_feasibility=true
        isaf = JuMP.is_solved_and_feasible(sub(model; index=index))
        sol_obj = jump_objective_all(sub(model; index=index); require_feasibility=false)

        model.info[:results][:subs][index][:obj] = isaf ? sol_obj.primal : missing
        model.info[:results][:subs][index][:obj_dual] = sol_obj.dual
        model.info[:results][:subs][index][:obj_lb] = isaf ? minimum(skipmissing(sol_obj)) : missing
        model.info[:results][:subs][index][:obj_ub] = isaf ? maximum(skipmissing(sol_obj)) : missing

        nothing
    end

    return nothing
end

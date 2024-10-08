function execute!(model::Benders.DecomposedModel, query::Benders.Query.ExtractResultsSub)
    tol = sqrt(eps(Float64))    # TODO

    @timeit model.timer "extract solution" begin
        jm = Benders.sub(model; index=query.index)
        res = model.info[:results][:subs][query.index]

        if has_attribute_type(model, Benders.CutTypeMISFSZ)
            # TODO: this still fails to properly estimate an UB ...

            π_0 = JuMP.value(jm.ext[:dualization_π_0])

            if π_0 <= -tol
                # The last cut was an "optimality" cut, so it is "feasible".
                # NOTE: Using `abs(...)` since `π_0` may have a negative sign (if coming from conic duality).
                res[:obj_val_primal] = res[:obj_val_dual] = (
                    JuMP.value(jm[:obj_base]) + 
                    abs(JuMP.value(jm[:obj_param]) / JuMP.value(jm.ext[:dualization_π_0]))
                )

                # TODO: NOTE/WRITING about the comment below

                # See ""
                #       http://www.dei.unipd.it/~fisch/papers/Benders_mis_extended_draft.pdf
                # at the end of Section 4 (page 10) for a note why this is a valid approach to calculate the original
                # dual's objective value.
            elseif sqrt(JuMP.value(jm[:obj_base])^2 + JuMP.value(jm[:obj_param])^2) <= tol
                # All dual variables are zero, so the modified dual is not unbounded anymore, indicating that the
                # original primal is feasible. This allows picking the currently set value of `θ` as valid (upper bound)
                # of the original primal's objective value.
                terms = values(jm[:obj_param_π_0].terms)
                res[:obj_val_primal] = res[:obj_val_dual] = isempty(terms) ? 0.0 : only(terms)
            else
                # The last cut was a "feasibility" cut, so the original dual is unbounded.
                res[:obj_val_primal] = missing
                res[:obj_val_dual] = +Inf
            end
        else
            # TODO: properly refactor this, like above, and then account for it everywhere

            # TODO: require_feasibility=true
            isaf = JuMP.is_solved_and_feasible(Benders.sub(model; index=query.index))
            sol_obj = jump_objective_all(Benders.sub(model; index=query.index); require_feasibility=false)

            model.info[:results][:subs][query.index][:obj] = isaf ? sol_obj.primal : missing
            model.info[:results][:subs][query.index][:obj_dual] = sol_obj.dual
            model.info[:results][:subs][query.index][:obj_lb] = isaf ? minimum(skipmissing(sol_obj)) : missing
            model.info[:results][:subs][query.index][:obj_ub] = isaf ? maximum(skipmissing(sol_obj)) : missing
        end
        nothing
    end

    return true
end

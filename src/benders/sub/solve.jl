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

    jm = Benders.sub(model; index = query.index)
    mis_fsz_cuts = has_attribute_type(model, Benders.CutTypeMISFSZ)

    # Fix the current main-model solution in the sub-model.
    @timeit model.timer "fix variables" begin
        if get(jm.ext, :dualization_is_dualized, false)
            # Dualized model.
            v2v_map = jm.ext[:dualization_var_to_vi]

            # Construct all objective parts.
            jm[:obj_base] = jm.ext[:dualization_obj_base]

            jm[:obj_param] = JuMP.AffExpr(0.0)
            jm[:obj_param_π_0] = JuMP.AffExpr(0.0)
            for elem in jm.ext[:dualization_obj_param]
                if haskey(v2v_map, elem[1])
                    JuMP.add_to_expression!(
                        jm[:obj_param],
                        model.info[:results][:main][:sol][v2v_map[elem[1]]] * elem[3],
                        elem[2],
                    )
                elseif mis_fsz_cuts
                    # The "artifical" π_0 variable.
                    JuMP.add_to_expression!(
                        jm[:obj_param_π_0],
                        JuMP.value(Benders.main(model)[:θ][query.index]) * elem[3],
                        jm.ext[:dualization_π_0],
                    )
                else
                    @error "Missing variable mapping for dualized variable"
                end
            end

            jm[:obj] = jm[:obj_param] + jm[:obj_base] + jm[:obj_param_π_0]

            # if mis_fsz_cuts
            #     # Update & add the MISFSZ objective term.
            #     θ_current = JuMP.value(Benders.main(model)[:θ][query.index])
            #     jm[:obj_misfsz] = jm[:π_0] * θ_current
            #     jm[:obj] -= jm[:obj_misfsz]
            # end

            # Update objective.
            JuMP.@objective(jm, Max, jm[:obj])

            # Base.Main.@infiltrate
        else
            # Standard (primal) model.
            for vi in jm[:y].axes[1]
                JuMP.fix(jm[:y][vi], model.info[:results][:main][:sol][vi]; force = true)
            end
        end
    end

    # Solve the sub's JuMP model.
    @timeit model.timer "optimize" JuMP.optimize!(jm)

    return nothing
end

# using JuMP
# import Gurobi

# JuMP.write_to_file(jm, "tmp.lp")
# jm = JuMP.read_from_file("tmp.mps")

# JuMP.set_optimizer(jm, Gurobi.Optimizer)

# JuMP.set_attribute(jm, "Presolve", 0)
# JuMP.set_attribute(jm, "Method", 2)
# JuMP.set_attribute(jm, "PreDual", 0)
# JuMP.set_attribute(jm, "NumericFocus", 3)
# JuMP.set_attribute(jm, "InfUnbdInfo", 1)
# JuMP.set_attribute(jm, "DualReductions", 0)

# JuMP.optimize!(jm)

# JuMP.solution_summary(jm)

# JuMP.constraint_by_name(jm, "cglp_normalization")

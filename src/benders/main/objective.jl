# TODO: transform the following functions into queries!

function get_obj_lb_main(model::Benders.DecomposedModel)
    # TODO: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592
    return jump_objective_lb(Benders.main(model))
end

function get_obj_ub_main(model::Benders.DecomposedModel)
    # TODO: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592
    return jump_objective_ub(Benders.main(model))
end

function get_obj_base_main(model::Benders.DecomposedModel)
    if !haskey(Benders.main(model), :obj_base)
        JuMP.@expression(
            Benders.main(model),
            obj_base,
            model.lpmd.c_offset + model.lpmd.c[model.vis[1]]' * Benders.main(model)[:x].data
        )
    end

    return Benders.main(model)[:obj_base]
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveConstOnly)
    @error "Currently not implemented"
    return false

    # m = Benders.main(model)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveDefault)
    JuMP.@expression(Benders.main(model), obj_full, get_obj_base_main(model) + sum(Benders.main(model)[:θ]))
    JuMP.@objective(Benders.main(model), Min, obj_full)

    return true
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveInCuts)
    @error "Currently not implemented"
    return false

    # m = bd_main(model)
    # idx_v = model.idx_model_vars[1]
    # @expression(m, obj, model.lpmd.c[idx_v]' * m[:x].data)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

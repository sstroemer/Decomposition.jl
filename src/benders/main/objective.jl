function get_obj_lb_main(model::Benders.DecomposedModel)
    # TODO: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592
    return jump_objective_lb(main(model))
end

function get_obj_ub_main(model::Benders.DecomposedModel)
    # TODO: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592
    return jump_objective_ub(main(model))
end

function get_obj_base_main(model::Benders.DecomposedModel)
    if !haskey(main(model), :obj_base)
        JuMP.@expression(main(model), obj_base, model.lpmd.c_offset + model.lpmd.c[model.vis[1]]' * main(model)[:x].data)
    end

    return main(model)[:obj_base]
end

function modify(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveConstOnly)
    @error "Currently not implemented"
    return nothing

    # m = main(model)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

function modify(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveDefault)
    JuMP.@expression(main(model), obj_full, get_obj_base_main(model) + sum(main(model)[:θ]))
    JuMP.@objective(main(model), Min, obj_full)

    add_attribute!(model, attribute)
    return nothing
end

function modify(model::Benders.DecomposedModel, attribute::Benders.Main.ObjectiveInCuts)
    @error "Currently not implemented"
    return nothing

    # m = bd_main(model)
    # idx_v = model.idx_model_vars[1]
    # @expression(m, obj, model.lpmd.c[idx_v]' * m[:x].data)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

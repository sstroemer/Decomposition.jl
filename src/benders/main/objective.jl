abstract type AbstractObjectiveType <: AbstractDecompositionAttribute end

@kwdef struct ObjectiveConstOnly <: AbstractObjectiveType; end
@kwdef struct ObjectiveDefault <: AbstractObjectiveType; end
@kwdef struct ObjectiveInCuts <: AbstractObjectiveType; end

function modify(model::DecomposedModel, attribute::ObjectiveConstOnly)
    @error "Currently not implemented"
    return nothing

    # m = main(model)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

function modify(model::DecomposedModel, attribute::ObjectiveDefault)
    JuMP.@objective(
        main(model),
        Min,
        model.lpmd.c_offset + model.lpmd.c[model.vis[1]]' * main(model)[:x].data + sum(main(model)[:θ])
    )

    add_attribute!(model.attributes, attribute)
    return nothing
end

function modify(model::DecomposedModel, attribute::ObjectiveInCuts)
    @error "Currently not implemented"
    return nothing

    # m = bd_main(model)
    # idx_v = model.idx_model_vars[1]
    # @expression(m, obj, model.lpmd.c[idx_v]' * m[:x].data)
    # @objective(m, Min, model.lpmd.c_offset + sum(m[:θ]))
    # push!(model.attributes, attribute)
    # return nothing
end

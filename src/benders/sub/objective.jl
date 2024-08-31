function modify(model::Benders.DecomposedModel, attribute::Benders.Sub.ObjectiveSelf)
    vis_main = model.vis[1]

    for i in 1:(length(model.models) - 1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = sub(model; index=i)
        vis_sub = model.vis[1 + i]
    
        m_sub[:obj] = JuMP.AffExpr(0.0)
        for vi in vis_sub
            (vi in vis_main) && continue
            JuMP.add_to_expression!(m_sub[:obj], model.lpmd.c[vi], m_sub[:x][vi])
        end

        JuMP.@objective(m_sub, Min, m_sub[:obj])
    end

    add_attribute!(model, attribute)
    return nothing
end

function modify(model::Benders.DecomposedModel, attribute::Benders.Sub.ObjectiveShared)
    # TODO: Implement this
    @error "This does not work yet (multiple subs need to correctly SHARE main-variable costs)!"

    # m = bd_sub(model; index=attribute.index)
    # idx_v = model.idx_model_vars[1 + attribute.index]   # TODO: account for "share of main variables, with multiple sub-models"

    # obj = AffExpr(0.0)
    # for i in idx_v
    #     # (i in model.idx_model_vars[1]) && continue    # TODO: divide here!
    #     add_to_expression!(obj, model.lpmd.c[i], m[:x][i])
    # end
    # @objective(m, Min, obj)

    # push!(model.attributes, attribute)

    return nothing
end

function modify(model::Benders.DecomposedModel, attribute::Benders.Main.VirtualSoftBounds)
    for vi in model.vis[1]
        if isfinite(attribute.lower) && !JuMP.has_lower_bound(main(model)[:x][vi])
            JuMP.set_lower_bound(main(model)[:x][vi], attribute.lower)
        end

        if isfinite(attribute.upper) && !JuMP.has_upper_bound(main(model)[:x][vi])
            JuMP.set_upper_bound(main(model)[:x][vi], attribute.upper)
        end
    end

    add_attribute!(model, attribute)
    return nothing
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Main.VirtualSoftBounds)
    jm = Benders.main(model)

    for vi in model.vis[1]
        if isfinite(attribute.lower) && !JuMP.has_lower_bound(jm[:x][vi])
            JuMP.set_lower_bound(jm[:x][vi], attribute.lower)
        end

        if isfinite(attribute.upper) && !JuMP.has_upper_bound(jm[:x][vi])
            JuMP.set_upper_bound(jm[:x][vi], attribute.upper)
        end
    end

    return true
end

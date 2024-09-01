function modify(model::Benders.DecomposedModel, attribute::Benders.Main.RegularizationLevelSet)
    if !has_attribute_type(model, Benders.Main.ObjectiveDefault)
        # TODO: fix this by implementing it for other ways
        @error "Missing ObjectiveDefault attribute"
        return nothing
    end

    # Add the "parameter variable" `L(α)`.
    JuMP.@variable(main(model), reg_levelset_L)

    # Add the level-set constraint.
    JuMP.@constraint(main(model), reg_levelset_constraint, main(model)[:obj_full] <= reg_levelset_L)

    # TODO: set UB = 1eX for iterations until a first upper bound is found

    add_attribute!(model, attribute)
    return nothing
end

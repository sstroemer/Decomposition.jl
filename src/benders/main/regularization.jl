function apply!(model::Benders.DecomposedModel, attribute::Benders.Main.RegularizationLevelSet)
    if !has_attribute_type(model, Benders.Main.ObjectiveDefault)
        # TODO: fix this by implementing it for other ways
        @error "Missing `ObjectiveDefault` attribute"
        return false
    end

    # Add the "parameter variable" `L(α)`.
    JuMP.@variable(Benders.main(model), reg_levelset_L)

    # Add the level-set constraint.
    JuMP.@constraint(Benders.main(model), reg_levelset_constraint, Benders.main(model)[:obj_full] <= reg_levelset_L)

    # TODO: set UB = 1eX for iterations until a first upper bound is found

    return true
end

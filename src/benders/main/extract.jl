function extract_main(model::Benders.DecomposedModel)
    if has_attribute_type(model, Benders.Main.RegularizationLevelSet)
        model.info[:results][:main][:obj_base] = JuMP.value(get_obj_base_main(model))
        # `model.info[:results][:main][:obj_lb]` was already set in solve
        model.info[:results][:main][:obj_ub] = get_obj_ub_main(model)
        model.info[:results][:main][:sol] = JuMP.value.(main(model)[:x])
    else
        model.info[:results][:main][:obj_base] = JuMP.value(get_obj_base_main(model))
        model.info[:results][:main][:obj_lb] = get_obj_lb_main(model) 
        model.info[:results][:main][:obj_ub] = get_obj_ub_main(model)
        model.info[:results][:main][:sol] = JuMP.value.(main(model)[:x])
    end

    return nothing
end

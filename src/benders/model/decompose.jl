function apply_annotations!(model::Benders.DecomposedModel)
    @info "Begin applying annotations; building models"
    # TODO: identify "trivial" sub-models here, and "remove" them

    # Prepare & save the full variable / constraint indices.
    vis_main = sort(union([v for (k, v) in model.annotations[:variables] if startswith(string(k), "main")]...))
    cis_main = sort(union([c for (k, c) in model.annotations[:constraints] if startswith(string(k), "main")]...))
    push!(model.vis, vis_main)       
    push!(model.cis, cis_main)

    # Create & save main-model.
    push!(model.models, Benders.lpmd_to_jump(model, vis_main, cis_main; name="main", optimizer=model.f_opt_main()))

    # TODO / NOTE for WRITING:
    # Grouping sub-models does not really work without extracting distinct cuts afterwards. It may be seen as a different take on "single-cut" (instead of multicut).
    # It potentially "hides" the effect of Y on sub-model S1, because S2 "dominates" the necessary decision. A work around would be to create the merged model with
    # duplicate variable copies: Instead of doing "unique", the variable "x_15111" should be created for both S1 and S2, their objective functions should just be added
    # together, and then the cuts can be separated as post-processing. This would allow to see the effect of Y on S1 and S2 separately, but still have the "grouped" effect.

    # Prepare correct index lists for each sub-model.
    n_sub_models = count((k) -> startswith(string(k), "sub_"), keys(model.annotations[:variables]))
    vis_sub = Vector{Int64}[]
    cis_sub = Vector{Int64}[]

    if !has_attribute_type(model, Benders.Sub.PreGroupModels)
        for i in 1:n_sub_models
            push!(vis_sub, model.annotations[:variables][Symbol("sub_$i")])
            push!(cis_sub, model.annotations[:constraints][Symbol("sub_$i")])
        end
    else
        cfg_nof_groups = get_attribute(model, Benders.Sub.PreGroupModels, :nof_groups)

        # Estimate problem "size / complexity" for each sub-model.
        problem_size = [
            sqrt(
                length(model.annotations[:variables][Symbol("sub_$i")])^2 +
                length(model.annotations[:constraints][Symbol("sub_$i")])^2
            ) for i in 1:n_sub_models
        ]

        # Sort the sub-models by their complexity.
        sorted_indices = sortperm(problem_size; rev=true)

        # Calculate the number of groups, either:
        nof_groups = (
            if cfg_nof_groups > 0
                # Based on the number that the user chose manually.
                min(n_sub_models, cfg_nof_groups)
            else
                # Estimate automatically, by trying to create "balanced" group sizes.
                est = problem_size[sorted_indices[1:2]]
                for i in sorted_indices[3:end]
                    if (problem_size[i] + est[end]) <= est[end - 1]
                        est[end] += problem_size[i]
                    else
                        push!(est, problem_size[i])
                    end
                end
                length(est)
            end
        )

        groups = [Int[] for _ in 1:nof_groups]
        group_size = zeros(Float64, nof_groups)
    
        # Create groups by greedily adding the largest sub-models first, and always to the group with the currently
        # smallest total size.
        for idx in sorted_indices
            smallest_group = argmin(group_size)
            push!(groups[smallest_group], idx)
            group_size[smallest_group] += problem_size[idx]
        end
    
        for group in groups
            # Multiple models may try to add the same variable copies (from main) to the same group.
            push!(vis_sub, unique([vi for g in group for vi in model.annotations[:variables][Symbol("sub_$g")]]))
            push!(cis_sub, unique([ci for g in group for ci in model.annotations[:constraints][Symbol("sub_$g")]]))
        end
    end

    for i in eachindex(vis_sub)
        # Save the full variable / constraint indices.
        push!(model.vis, vis_sub[i])
        push!(model.cis, cis_sub[i])

        # Create & save sub-model.
        push!(model.models, Benders.lpmd_to_jump(model, vis_sub[i], cis_sub[i]; name="sub_$i", optimizer=model.f_opt_sub()))
    end

    # Cache the variable indices in main that link to a specific sub-model.
    set_vis_main = Set(vis_main)
    model.annotations[:linking_variables] = Dict()
    for i in eachindex(vis_sub)
        # TODO: relax the dict typing in the model definition to directly use `i` here as key
        model.annotations[:linking_variables][Symbol(i)] = intersect(Set(vis_sub[i]), set_vis_main)
    end

    # Create empty results dictionary.
    model.info[:results] = OrderedDict{Symbol, Any}(
        :main => OrderedDict{Symbol, Any}(),
        :subs => [OrderedDict{Symbol, Any}() for _ in Benders.subs(model)],
    )

    @info "Annotations successfully applied" n_main_models = 1 n_sub_models = length(vis_sub)
    return true
end

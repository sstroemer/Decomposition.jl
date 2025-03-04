function _preprocess_cuts_drop_non_binding!(model::Benders.DecomposedModel)
    !has_attribute_type(model, Benders.CutPostprocessingDropNonBinding) && return nothing

    @timeit model.timer "drop non-binding" begin
        attr = get_attribute(model, Benders.CutPostprocessingDropNonBinding)

        to_delete = []

        @timeit model.timer "scan" begin
            for k in keys(model.cuts)
                for cut in model.cuts[k]
                    if !haskey(cut.stats, :nonbinding)
                        cut.stats[:nonbinding] = 0
                    end

                    # Skip cuts that were already removed.
                    (cut.stats[:nonbinding] == -1) && continue

                    # Update the non-binding counter.
                    if JuMP.value(cut.cut_con) < attr.threshold
                        cut.stats[:nonbinding] = 0
                    else
                        cut.stats[:nonbinding] += 1
                    end

                    # Remove the cut if it was non-binding for a certain number of iterations.
                    if cut.stats[:nonbinding] >= attr.iterations
                        cut.stats[:nonbinding] = -1
                        push!(to_delete, cut)
                    end
                end
            end
        end

        # Remove the non-binding cuts.
        @timeit model.timer "delete" begin
            for cut in to_delete
                JuMP.delete(JuMP.owner_model(cut.cut_con), cut.cut_con)
            end
        end
    end
    
    return nothing
end

function postprocess_cuts!(model::Benders.DecomposedModel)
    _preprocess_cuts_drop_non_binding!(model)

    return nothing
end

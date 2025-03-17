function _binding_constraints(cuts; τ::Float64)
    binding_cons = JuMP.ConstraintRef[]

    for con in cuts
        (abs(JuMP.dual(con)) < τ) || continue
        push!(binding_cons, con)
    end

    return binding_cons
end

function _preprocess_cuts_drop_non_binding!(model::Benders.DecomposedModel)
    !has_attribute_type(model, Benders.CutPostprocessingDropNonBinding) && return nothing

    @timeit model.timer "drop non-binding" begin
        attr = get_attribute(model, Benders.CutPostprocessingDropNonBinding)

        to_delete = JuMP.JuMP.ConstraintRef{
            JuMP.Model,
            JuMP.MOI.ConstraintIndex{JuMP.MOI.ScalarAffineFunction{Float64}, JuMP.MOI.GreaterThan{Float64}},
            JuMP.ScalarShape,
        }[]

        @timeit model.timer "scan" begin
            for cuts in values(model.cuts)
                for cut in cuts::Vector{Decomposition.Benders.AbstractGeneralCut}
                    cs = cut.stats::Dict{Symbol, Any}
                    cnb = cs[:nonbinding]::Int64

                    # Skip cuts that were already removed.
                    (cnb < 0) && continue

                    con = cut.cut_con::JuMP.ConstraintRef{
                        JuMP.Model,
                        JuMP.MOI.ConstraintIndex{JuMP.MOI.ScalarAffineFunction{Float64}, JuMP.MOI.GreaterThan{Float64}},
                        JuMP.ScalarShape,
                    }

                    # Update the non-binding counter.
                    cs[:nonbinding] = (abs(JuMP.dual(con)) < attr.threshold ? (cnb + 1) : 0)::Int64

                    # Remove the cut if it was non-binding for a certain number of iterations.
                    if cs[:nonbinding]::Int64 >= attr.iterations
                        cs[:nonbinding] = (-1)::Int64
                        push!(to_delete, con)
                    end
                end
            end
        end

        # Remove the non-binding cuts.
        @timeit model.timer "delete" begin
            isempty(to_delete) || JuMP.delete(Benders.main(model)::JuMP.Model, to_delete)
        end
    end

    return nothing
end

function postprocess_cuts!(model::Benders.DecomposedModel)
    _preprocess_cuts_drop_non_binding!(model)

    return nothing
end

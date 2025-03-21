function _preprocess_cuts_remove_redundant!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    !has_attribute_type(model, Benders.CutPreprocessingRemoveRedundant) && return nothing

    attr = get_attribute(model, Benders.CutPreprocessingRemoveRedundant)
    tol = (rel_ctol = attr.rtol_coeff, rel_btol = attr.rtol_const)

    # Note: Optimality cuts contain θ and can therefore only be compared against once from the same sub-model.

    # TODO: why does this seem to hurt solver performance?

    total_valid_cuts = 0
    total_raw_cuts = 0

    opt_sub_map = Dict(i => [] for i in 1:(length(model.models)-1))
    for c in reverse(model.cuts[:optimality])
        push!(opt_sub_map[c.sub_model_index], c.cut_exp)
    end

    for cut_type in [:feasibility, :optimality, :misfsz]
        valid_cuts = []
        for cut in new_cuts[cut_type]
            if cut_type == :optimality
                any(jump_expressions_equal(cut[2], other; tol...) for other in opt_sub_map[cut[1]]) && continue
            else
                any(jump_expressions_equal(cut[2], other.cut_exp; tol...) for other in model.cuts[cut_type]) && continue
            end
            push!(valid_cuts, cut)
        end

        total_valid_cuts += length(valid_cuts)
        total_raw_cuts += length(new_cuts[cut_type])

        new_cuts[cut_type] = valid_cuts
    end

    # TODO: make that configurable
    τ = 0.25  # reduce by how much?
    π = 0.50  # when?
    if total_valid_cuts < total_raw_cuts * π
        for (i, dac) in enumerate(model.attributes)
            if dac.attribute == attr
                model.attributes[i] = DecompositionAttributeContainer(;
                    model = dac.model,
                    attribute = Benders.CutPreprocessingRemoveRedundant(;
                        rtol_coeff = tol.rel_ctol * τ,
                        rtol_const = tol.rel_btol * τ,
                    ),
                    iteration = dac.iteration,
                    dirty = dac.dirty,
                )
                break
            end
        end
    end

    return nothing
end

function _preprocess_cuts_stabilize_numerical_range!(
    model::Benders.DecomposedModel,
    new_cuts::Dict{Symbol, Vector{Any}},
)
    !has_attribute_type(model, Benders.CutPreprocessingStabilizeNumericalRange) && return nothing
    attr = get_attribute(model, Benders.CutPreprocessingStabilizeNumericalRange)

    # Remove extremely low coefficients.
    # Example: (4, -0.0016145468734458735 x[121] - 765000.7505230922 x[363] ... + 2.5318715674811035e10)

    for cut in new_cuts[:optimality]
        constant = cut[2].constant
        constant_factor = abs(constant) ./ abs.(cut[2].terms.vals)

        for (var, coeff, fact) in zip(cut[2].terms.keys, cut[2].terms.vals, constant_factor)
            if fact > attr.const_factor_threshold
                max_delta = Inf

                if coeff > 0 && JuMP.has_lower_bound(var)
                    max_delta = JuMP.lower_bound(var) * coeff
                elseif coeff < 0 && JuMP.has_upper_bound(var)
                    max_delta = JuMP.upper_bound(var) * coeff
                end

                if abs(max_delta / constant) < attr.const_factor_elimination_max_rel_delta
                    delete!(cut[2].terms, var)
                    cut[2].constant += max_delta
                end
            end
        end
    end

    return nothing
end

function _preprocess_cuts_make_unique!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    !has_attribute_type(model, Benders.CutPreprocessingMakeUnique) && return nothing

    attr = get_attribute(model, Benders.CutPreprocessingMakeUnique)
    tol = (rel_ctol = attr.rtol_coeff, rel_btol = attr.rtol_const)

    # Note: Optimality cuts contain θ and can therefore only be compared against once from the same sub-model.

    for cut_type in [:feasibility, :optimality, :misfsz]
        valid_cuts = []
        for cut in new_cuts[cut_type]
            if cut_type == :optimality
                any(jump_expressions_equal(cut[2], other[2]; tol...) for other in valid_cuts if cut[1] == other[1]) &&
                    continue
            else
                any(jump_expressions_equal(cut[2], other[2]; tol...) for other in valid_cuts) && continue
            end
            push!(valid_cuts, cut)
        end
        new_cuts[cut_type] = valid_cuts
    end

    return nothing
end

function preprocess_cuts!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    # TODO: keep stats on number of preprocessing steps, etc.
    _preprocess_cuts_make_unique!(model, new_cuts)
    _preprocess_cuts_remove_redundant!(model, new_cuts)
    # _preprocess_cuts_stabilize_numerical_range!(model, new_cuts)

    return nothing
end

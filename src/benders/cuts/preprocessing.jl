function _preprocess_cuts_remove_redundant!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    !has_attribute_type(model, Benders.CutPreprocessingRemoveRedundant) && return nothing

    # TODO: feasibility cuts can also be redundant within an iteration/ between sub-models: (1, -x[2905] + 31965.28), (4, -x[2905] + 36838.016), ...
    for type in [:feasibility, :optimality]
        valid = [true for _ in new_cuts[type]]

        for i in eachindex(new_cuts[type])
            cut = new_cuts[type][i]
            
            iter_oldcuts = (c for c in model.cuts[type] if c.sub_model_index == cut[1])
            isempty(iter_oldcuts) && continue

            last_opt_cut = last(c for c in model.cuts[type] if c.sub_model_index == cut[1])
            cut[2] != last_opt_cut.cut_exp && continue

            # This cut is identical to an old one.
            valid[i] = false
        end

        new_cuts[type] = [new_cuts[type][i] for i in eachindex(new_cuts[type]) if valid[i]]
    end

    return nothing
end

function _preprocess_cuts_stabilize_numerical_range!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
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

function preprocess_cuts!(model::Benders.DecomposedModel, new_cuts::Dict{Symbol, Vector{Any}})
    # TODO: keep stats on number of preprocessing steps, etc.
    _preprocess_cuts_remove_redundant!(model, new_cuts)
    _preprocess_cuts_stabilize_numerical_range!(model, new_cuts)

    return nothing
end

abstract type AbstractCut <: Decomposition_AbtractCut end

@kwdef struct GeneralOptimalityCut <: AbstractCut
    iteration::Int
    sub_model_index::Int = -1
    
    cut_exp::JuMP.AffExpr
    cut_con::JuMP.ConstraintRef

    stats::Dict = Dict(
        :abs_constant => abs(cut_exp.constant),
        :range_coeff => extrema(abs.(cut_exp.terms.vals); init=(0.0, 0.0))
    )
end

@kwdef struct GeneralFeasibilityCut <: AbstractCut
    iteration::Int
    sub_model_index::Int = -1
    
    cut_exp::JuMP.AffExpr
    cut_con::JuMP.ConstraintRef

    stats::Dict = Dict(
        :abs_constant => abs(cut_exp.constant),
        :range_coeff => extrema(abs.(cut_exp.terms.vals); init=(0.0, 0.0))
    )
end

# =========================================================================
# Cut processing

abstract type AbstractCutProcessing <: AbstractDecompositionAttribute end

@kwdef struct CutPreprocessingRemoveRedundant <: AbstractCutProcessing end

@kwdef struct CutPreprocessingStabilizeNumericalRange <: AbstractCutProcessing
    const_factor_threshold::Float64                  # e.g., 1e10
    const_factor_elimination_max_rel_delta::Float64  # e.g., 1e-4
end

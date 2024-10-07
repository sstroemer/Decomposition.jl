# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: General Abstract Cut Types                                                                             │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

abstract type AbstractCutType <: AbstractGeneralAttribute end
abstract type AbstractFeasibilityCutType <: AbstractCutType end
abstract type AbstractOptimalityCutType <: AbstractCutType end

abstract type AbstractCutProcessing <: AbstractGeneralAttribute end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Cut Types                                                                                              │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct FeasibilityCutTypeSingle <: AbstractFeasibilityCutType end
@kwdef struct FeasibilityCutTypeMulti <: AbstractFeasibilityCutType end
@kwdef struct FeasibilityCutTypeAggregated <: AbstractFeasibilityCutType end
@kwdef struct FeasibilityCutTypeAdaptive <: AbstractFeasibilityCutType end

@kwdef struct OptimalityCutTypeSingle <: AbstractOptimalityCutType end
@kwdef struct OptimalityCutTypeMulti <: AbstractOptimalityCutType end
@kwdef struct OptimalityCutTypeAggregated <: AbstractOptimalityCutType end
@kwdef struct OptimalityCutTypeAdaptive <: AbstractOptimalityCutType end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ CUTS: Optimality                                                                                                   │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct GeneralOptimalityCut <: AbstractGeneralCut
    iteration::Int
    sub_model_index::Int = -1
    
    cut_exp::JuMP.AffExpr
    cut_con::JuMP.ConstraintRef

    stats::Dict = Dict(
        :abs_constant => abs(cut_exp.constant),
        :range_coeff => extrema(abs.(cut_exp.terms.vals); init=(0.0, 0.0))
    )
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ CUTS: Feasibility                                                                                                  │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct GeneralFeasibilityCut <: AbstractGeneralCut
    iteration::Int
    sub_model_index::Int = -1
    
    cut_exp::JuMP.AffExpr
    cut_con::JuMP.ConstraintRef

    stats::Dict = Dict(
        :abs_constant => abs(cut_exp.constant),
        :range_coeff => extrema(abs.(cut_exp.terms.vals); init=(0.0, 0.0))
    )
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Cut Processing                                                                                         │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct CutPreprocessingRemoveRedundant <: AbstractCutProcessing end

@kwdef struct CutPreprocessingStabilizeNumericalRange <: AbstractCutProcessing
    const_factor_threshold::Float64                  # e.g., 1e10
    const_factor_elimination_max_rel_delta::Float64  # e.g., 1e-4
end

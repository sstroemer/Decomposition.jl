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

"""
    CutTypeMISFSZ

"Minimal Infeasible Subsystem" (MIS) cuts, based on the paper by Fischetti, Salvagnin, & Zanette (FSZ).
DOI 10.1007/s10107-010-0365-7
"""
@kwdef struct CutTypeMISFSZ <: AbstractCutType end

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

# TODO: add the comparison tolerances here as parameters
@kwdef struct CutPreprocessingRemoveRedundant <: AbstractCutProcessing
    rtol_coeff::Float64 = eps(Float64)
    rtol_const::Float64 = eps(Float64)
end

@kwdef struct CutPreprocessingMakeUnique <: AbstractCutProcessing
    rtol_coeff::Float64 = eps(Float64)
    rtol_const::Float64 = eps(Float64)
end

@kwdef struct CutPreprocessingStabilizeNumericalRange <: AbstractCutProcessing
    const_factor_threshold::Float64                  # e.g., 1e10
    const_factor_elimination_max_rel_delta::Float64  # e.g., 1e-4
end

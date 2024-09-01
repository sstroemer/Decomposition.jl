module Main

import ...Decomposition: AbstractDecompositionAttribute

abstract type AbstractCutType <: AbstractDecompositionAttribute end

abstract type AbstractFeasibilityCutType <: AbstractCutType end
abstract type AbstractOptimalityCutType <: AbstractCutType end

abstract type AbstractRegularization <: AbstractDecompositionAttribute; end

@kwdef struct FeasibilityCutTypeSingle <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeMulti <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeAggregated <: AbstractFeasibilityCutType; end
@kwdef struct FeasibilityCutTypeAdaptive <: AbstractFeasibilityCutType; end

@kwdef struct OptimalityCutTypeSingle <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeMulti <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeAggregated <: AbstractOptimalityCutType; end
@kwdef struct OptimalityCutTypeAdaptive <: AbstractOptimalityCutType; end

abstract type AbstractObjectiveType <: AbstractDecompositionAttribute end

@kwdef struct ObjectiveConstOnly <: AbstractObjectiveType; end
@kwdef struct ObjectiveDefault <: AbstractObjectiveType; end
@kwdef struct ObjectiveInCuts <: AbstractObjectiveType; end

@kwdef struct RegularizationLevelSet <: AbstractRegularization
    alpha::Float64
    infeasible_alpha_step::Float64

    safety_max_infeasible_resolve::Int = 10
    
    _internal::Dict = Dict() # TODO: convert to proper fields
end

@kwdef struct RegularizationTrustRegion <: AbstractRegularization
    norm::Symbol # l1, l2, inf, ...
end

@kwdef struct VirtualSoftBounds <: AbstractDecompositionAttribute
    lower::Float64 = -Inf
    upper::Float64 = +Inf
end

end

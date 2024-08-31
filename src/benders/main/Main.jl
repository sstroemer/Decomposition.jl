module Main

import ...Decomposition: AbstractDecompositionAttribute

abstract type AbstractCutType <: AbstractDecompositionAttribute end

abstract type AbstractFeasibilityCutType <: AbstractCutType end
abstract type AbstractOptimalityCutType <: AbstractCutType end

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

@kwdef struct VirtualSoftBounds <: AbstractDecompositionAttribute
    lower::Float64 = -Inf
    upper::Float64 = +Inf
end

end

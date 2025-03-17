module Main

import ..AbstractGeneralAttribute as _AGA

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: General Abstract Types                                                                                 │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

abstract type AbstractMainModelAttribute <: _AGA end

abstract type AbstractRegularization <: AbstractMainModelAttribute end
abstract type AbstractObjectiveType <: AbstractMainModelAttribute end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Objective Function                                                                                     │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct ObjectiveConstOnly <: AbstractObjectiveType end
@kwdef struct ObjectiveDefault <: AbstractObjectiveType end
@kwdef struct ObjectiveInCuts <: AbstractObjectiveType end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Regularization                                                                                         │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct RegularizationLevelSet <: AbstractRegularization
    alpha::Float64
    infeasible_alpha_step::Float64

    safety_max_infeasible_resolve::Int = 10

    _internal::Dict = Dict() # TODO: convert to proper fields
end

@kwdef struct RegularizationTrustRegion <: AbstractRegularization
    norm::Symbol # l1, l2, inf, ...
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Auxiliary                                                                                              │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct VirtualSoftBounds <: AbstractMainModelAttribute
    lower::Float64 = -Inf
    upper::Float64 = +Inf
end

end

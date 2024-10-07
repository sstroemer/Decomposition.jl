module Sub

import JuMP
import ..AbstractGeneralAttribute as _AGA

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: General Abstract Types                                                                                 │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

abstract type AbstractSubModelAttribute <: _AGA end
abstract type AbstractFeasibilityMode <: AbstractSubModelAttribute end
abstract type AbstractObjectiveType <: AbstractSubModelAttribute end
abstract type AbstractRelaxationType <: AbstractSubModelAttribute end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Objective Function                                                                                     │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct ObjectiveSelf <: AbstractObjectiveType
    index::Int64 = -1
end

@kwdef struct ObjectiveShared <: AbstractObjectiveType
    index::Int64 = -1
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Relaxation                                                                                             │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct RelaxationAll <: AbstractRelaxationType
    index::Int64 = -1
    penalty::Float64 = 1e7

    _penalty_map::Dict{JuMP.ConstraintRef, JuMP.AffExpr} = Dict{JuMP.ConstraintRef, JuMP.AffExpr}()
end

@kwdef struct RelaxationLinked <: AbstractRelaxationType
    index::Int64 = -1
    penalty::Float64 = 1e7
end

@kwdef struct RelaxationRegex <: AbstractRelaxationType
    index::Int64 = -1
    penalty::Float64 = 1e7

    regex::Regex
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Auxiliary                                                                                              │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct PreGroupModels <: AbstractSubModelAttribute
    nof_groups::Int64 = -1
end

end

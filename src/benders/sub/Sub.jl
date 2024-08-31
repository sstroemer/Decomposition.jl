module Sub

import JuMP
import ...Decomposition: AbstractDecompositionAttribute, AbstractDecompositionQuery

abstract type AbstractFeasibilityMode <: AbstractDecompositionAttribute end

abstract type AbstractObjectiveType <: AbstractDecompositionAttribute end

@kwdef struct ObjectiveSelf <: AbstractObjectiveType
    index::Int64 = -1
end

@kwdef struct ObjectiveShared <: AbstractObjectiveType
    index::Int64 = -1
end

abstract type AbstractRelaxationType <: AbstractDecompositionAttribute end

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

@kwdef struct QueryRelaxation <: AbstractDecompositionQuery; end

end

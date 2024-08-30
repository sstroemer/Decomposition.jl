module Sub

import JuMP

import ...Decomposition: AbstractDecompositionAttribute, AbstractDecompositionQuery
import ..Benders: DecomposedModel

abstract type AbstractSolverAttribute <: AbstractDecompositionAttribute end

abstract type AbstractFeasibilityMode <: AbstractDecompositionAttribute end

include("general.jl")
include("objective.jl")
include("relaxation.jl")
include("update.jl")

end

import .Sub

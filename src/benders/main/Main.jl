module Main

import JuMP

import ...Decomposition: AbstractDecompositionAttribute
import ..Benders: DecomposedModel

abstract type AbstractSolverAttribute <: AbstractDecompositionAttribute end

abstract type AbstractCutType <: AbstractDecompositionAttribute end

include("general.jl")
include("cuts.jl")

end

import .Main

module Solver

import JuMP

import ..Decomposition: AbstractDecompositionAttribute

abstract type AbstractSolverAttribute <: AbstractDecompositionAttribute end

@kwdef struct AlgorithmSimplex <: AbstractSolverAttribute
    model::Any
    mode::Symbol = :dual
end

@kwdef struct AlgorithmIPM <: AbstractSolverAttribute
    model::Any
    mode::Symbol = :none
    crossover::Bool = false
end

@kwdef struct ExtractDualRay <: AbstractSolverAttribute
    activate::Bool = true
end

include("gurobi.jl")
include("highs.jl")

function modify(jump_model::JuMP.Model, attribute::AbstractSolverAttribute)
    solver_functions = Dict(
        "Gurobi" => _gurobi_modify_jump,
        "HiGHS" => _highs_modify_jump,
    )

    solver = JuMP.solver_name(jump_model)
    
    if !haskey(solver_functions, solver)
        @error "Solver `$(solver)` is currently not supported"
    elseif !solver_functions[solver](jump_model, attribute)
        @error "Cannot set solver attribute" solver attribute
    end

    return nothing
end

end

import .Solver

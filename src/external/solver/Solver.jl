module Solver

import JuMP

import ..Decomposition: AbstractDecompositionAttribute

abstract type AbstractAttribute <: AbstractDecompositionAttribute end

@kwdef struct AlgorithmSimplex <: AbstractAttribute
    model::Any
    mode::Symbol = :dual
end

@kwdef struct AlgorithmIPM <: AbstractAttribute
    model::Any
    mode::Symbol = :none
    crossover::Bool = false
end

@kwdef struct ExtractDualRay <: AbstractAttribute
    model::Any
    activate::Bool = true
end

include("gurobi.jl")
include("highs.jl")

function _modify_jump(jump_model::JuMP.Model, attribute::AbstractAttribute)
    solver_functions = Dict(
        "Gurobi" => _gurobi_modify_jump,
        "HiGHS" => _highs_modify_jump,
    )

    solver = JuMP.solver_name(jump_model)
    
    if !haskey(solver_functions, solver)
        @error "Solver is currently not supported for this attribute" solver attribute
        return false
    elseif !solver_functions[solver](jump_model, attribute)
        @error "Cannot set solver attribute" solver attribute
        return false
    end

    return true
end

end

import .Solver

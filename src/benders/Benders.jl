module Benders

import JuMP
import Dualization
import OrderedCollections: OrderedDict

import ..Decomposition.LinearAlgebra: mat_vec_scalar
import ..Decomposition: TimerOutput, @timeit
import ..Decomposition: has_attribute_type, get_attribute, get_attributes, cache_get
import ..Decomposition: current_iteration
import ..Decomposition: jump_objective_lb, jump_objective_ub, jump_objective_all
import ..Decomposition: annotate!, apply!, execute!

import ..Decomposition as _ExtModDecomposition
import ..Decomposition.Solver as _ExtModSolver

const MOI = JuMP.MOI

abstract type AbstractGeneralAttribute <: _ExtModDecomposition.AbstractDecompositionAttribute end
abstract type AbstractGeneralCut <: _ExtModDecomposition.AbstractCut end
abstract type AbstractGeneralQuery <: _ExtModDecomposition.AbstractDecompositionQuery end

include("_definitions/model.jl")
include("model/lpmd_to_jump.jl")
include("_definitions/definitions.jl")

main(model::Benders.DecomposedModel) = model.models[1]
sub(model::Benders.DecomposedModel; index::Int) = model.models[1+index]
subs(model::Benders.DecomposedModel) = model.models[2:end]

function check_cut_type(jump_model::JuMP.Model; verbose::Bool = true)
    if JuMP.is_solved_and_feasible(jump_model)
        return :optimality
    end

    if JuMP.dual_status(jump_model) != JuMP.MOI.INFEASIBILITY_CERTIFICATE
        verbose && (@error "Turn off presolve, or any setting blocking extraction of dual rays")
        return :error
    end

    return :feasibility
end

end

import .Benders

include("model/model.jl")
include("model/decompose.jl")

include("main/main.jl")
include("sub/sub.jl")

include("cuts/cuts.jl")

include("util/util.jl")
include("iterate.jl")

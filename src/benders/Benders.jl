module Benders

import JuMP
import OrderedCollections: OrderedDict

import ..Decomposition.LinearAlgebra: mat_vec_scalar
import ..Decomposition: TimerOutput, @timeit
import ..Decomposition: has_attribute_type, get_attribute, get_attributes
import ..Decomposition: current_iteration
import ..Decomposition: jump_objective_lb, jump_objective_ub, jump_objective_all
import ..Decomposition: annotate!, apply!, execute!

import ..Decomposition as _DecompositionMainModule

const MOI = JuMP.MOI

abstract type AbstractGeneralAttribute <: _DecompositionMainModule.AbstractDecompositionAttribute end
abstract type AbstractGeneralCut <: _DecompositionMainModule.AbstractCut end
abstract type AbstractGeneralQuery <: _DecompositionMainModule.AbstractDecompositionQuery end

include("_definitions/model.jl")
include("model/model.jl")
include("_definitions/definitions.jl")

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

function apply!(model::Benders.DecomposedModel, attribute::Solver.AbstractAttribute)
    models = attribute.model == :main ? [Benders.main(model)] : Benders.subs(model)
    
    ret = true
    for m in models
        ret &= Solver._modify_jump(m, attribute)
    end
    
    return ret
end

include("main/main.jl")
include("sub/sub.jl")

include("cuts/cuts.jl")

include("util/util.jl")
include("iterate.jl")

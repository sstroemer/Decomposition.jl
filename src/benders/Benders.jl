module Benders

import JuMP
import OrderedCollections: OrderedDict

import ..Decomposition: AbstractDecompositionAttribute, AbstractDecomposedModel, TimerOutput

include("model.jl")

include("main/Main.jl")
include("sub/Sub.jl")
include("termination/Termination.jl")

import .Main
import .Sub
import .Termination

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

include("general.jl")

function has_attribute_type(model::Benders.DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    return any(attr isa type for attr in model.attributes)
end

function get_attributes(model::Benders.DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    attributes = AbstractDecompositionAttribute[]
    for attr in model.attributes
        attr isa type && push!(attributes, attr)
    end

    if isempty(attributes)
        @error "Model does not contain an attribute of type `$(T)`, check first with `has_attribute_type(...)`"
    end

    return attributes
end

function get_attribute(model::Benders.DecomposedModel, type::Type{T}) where T <: AbstractDecompositionAttribute
    attributes = get_attributes(model, type)

    if length(attributes) > 1
        @error "Attributes of type `$(T)` not unique, use `get_attributes(...)` instead"
        return nothing
    end

    return attributes[1]
end

function modify(model::Benders.DecomposedModel, attribute::Solver.AbstractAttribute)
    models = attribute.model == :main ? [main(model)] : subs(model)
    for m in models
        Solver._modify_jump(m, attribute)
    end
    add_attribute!(model, attribute)
    return nothing
end

include("main/cuts.jl")
include("main/general.jl")
include("main/objective.jl")
include("main/regularization.jl")
include("main/solve.jl")    # should this be inside Main? (Main.solve, instead of solve_main)
include("main/extract.jl")  # should this be inside Main? (Main.extract, instead of extract_main)

include("sub/objective.jl")
include("sub/relaxation.jl")

include("termination/termination.jl")
include("iterate.jl")

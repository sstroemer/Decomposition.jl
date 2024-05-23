abstract type AbstractDecompositionNode <: AbstractNode end

include("admm.jl")
include("benders.jl")

function decompose(node::AbstractNode; method::Symbol = :benders_simple)
    @debug "decompose(::AbstractNode; ::Symbol)" node method

    method == :benders_simple && return SimpleBendersDecomposition(node)

    @critical "Unknown decomposition method" method
end

decompose(node::String; method::Symbol = :benders_simple) = decompose(from_file(node); method = method)

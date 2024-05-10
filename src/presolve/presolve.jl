include("papilo.jl")

function presolve(node::AbstractNode; algorithm::Symbol = :papilo)
    @debug "presolve(::AbstractNode; ::Symbol)" node

    (algorithm == :papilo) && (return _presolve_papilo(node))

    @critical "`presolve` got unknown algorithm, currently supported: `:papilo`" algorithm
end

presolve(node::AbstractFileNode; algorithm::Symbol = :papilo) = presolve(to_model(node); algorithm = algorithm)
presolve(node::String; algorithm::Symbol = :papilo) = presolve(from_file(node); algorithm = algorithm)

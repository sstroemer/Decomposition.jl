include("papilo.jl")


function presolve(node::AbstractNode; algorithm::Symbol = :papilo)
    @debug "presolve(::AbstractNode; ::Symbol)" node
    
    (algorithm == :papilo) && (return _presolve_papilo(node))
    
    @critical "`presolve` got unknown algorithm, currently supported: `:papilo`" algorithm
end

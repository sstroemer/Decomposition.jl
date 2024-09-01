module Termination

import ...Decomposition: AbstractDecompositionAttribute

abstract type AbstractCriterion <: AbstractDecompositionAttribute end

@kwdef struct Stop <: AbstractCriterion
    iterations::Int64 = -1
    
    seconds::Real = -1

    opt_gap_abs::Real = -1
    opt_gap_rel::Real = -1
end

end

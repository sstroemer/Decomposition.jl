module Termination

import ..AbstractGeneralAttribute as _AGA

abstract type AbstractGeneralCriterion <: _AGA end

@kwdef struct Stop <: AbstractGeneralCriterion
    iterations::Int64 = -1
    
    seconds::Real = -1

    opt_gap_abs::Real = -1
    opt_gap_rel::Real = -1
end

end

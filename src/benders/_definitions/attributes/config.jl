module Config

import ..AbstractGeneralAttribute as _AGA

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: General Abstract Types                                                                                 │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

abstract type AbstractProblemConfiguration <: _AGA end
abstract type AbstractModelConfiguration <: _AGA end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Problem Configuration                                                                                  │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

"""
    TotalTimesteps(T::Int64)

Problem configuration for the total number of timesteps in the original problem.

# Arguments
- `T::Int64`: Total number of timesteps in the original problem.
"""
@kwdef struct TotalTimesteps <: AbstractProblemConfiguration
    T::Int64
end

"""
    NumberOfTemporalBlocks(n::Int64 = 1)

Problem configuration for the number of temporal blocks in the original problem.

# Arguments
- `n::Int64 = 1`: Number of temporal blocks in the original problem. The default `1` means that the whole problem is
  considered as one continuous block.
"""
@kwdef struct NumberOfTemporalBlocks <: AbstractProblemConfiguration
    n::Int64 = 1
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ ATTRIBUTES: Model Configuration                                                                                    │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

"""
    ModelVerbosity(verbosity::Int64 = 2)

Model configuration for the verbosity of the model.

# Arguments
- `verbosity::Int64 = 2`: Verbosity level of the model. Possible values (always including all of the lower levels):
    - `0`: silent
    - `1`: active logging >= warning
    - `2`: output iteration print log
    - `3`: active logging >= info
    - `4`: active logging >= debug
"""
@kwdef struct ModelVerbosity <: AbstractModelConfiguration
    verbosity::Int64 = 2
end

"""
    ModelDebug(enable::Bool = false)

Model configuration for the debug mode of the model.

# Arguments
- `debug::Bool = false`: Debug mode of the model.
"""
@kwdef struct ModelDebug <: AbstractModelConfiguration
    enable::Bool = false
end

"""
    ModelDirectMode(enable::Bool)

Model configuration for using `JuMP.direct_model` instead of `JuMP.Model`.

# Arguments
- `enable::Bool`: Direct mode of the model.
"""
@kwdef struct ModelDirectMode <: AbstractModelConfiguration
    enable::Bool = true
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ DEFAULTS                                                                                                           │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

union!(parentmodule(@__MODULE__)._ATTRIBUTE_DEFAULTS, [
    NumberOfTemporalBlocks(),
    ModelVerbosity(),
    ModelDebug(),
    ModelDirectMode(),
])

end

module Query

import ..AbstractGeneralQuery as _AGQ

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ QUERIES: General Abstract Type                                                                                     │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

abstract type AbstractExecuteQuery <: _AGQ end
abstract type AbstractResultQuery <: _AGQ end
abstract type AbstractQualityQuery <: _AGQ end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ QUERIES: Execute                                                                                                   │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct SolveMain <: AbstractExecuteQuery end

@kwdef struct SolveSub <: AbstractExecuteQuery
    index::Int64
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ QUERIES: Results                                                                                                   │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct ExtractResultsMain <: AbstractResultQuery end

@kwdef struct ExtractResultsSub <: AbstractResultQuery
    index::Int64
end

# ╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
# │ QUERIES: Quality                                                                                                   │
# ╰ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶ ╶╯

@kwdef struct Relaxation <: AbstractQualityQuery end

end

import .Query

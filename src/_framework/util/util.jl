import SparseArrays
using Suppressor

import JuMP
import Dualization
const MOI = JuMP.MOI

import SCIP_PaPILO_jll
import HiGHS    # TODO: remove this after allowing to set an optimizer in the root node

import Printf

macro critical(msg, args...)
    msg = string(msg)
    return esc(quote
        @error $msg $(args...)
        error($msg)
    end)
end

function _norm_strname(str::String, width::Int64 = 50)
    textwidth(str) <= width && return rpad(str, width)

    prefix_length = div(width - 3, 2) + div(width - 3, 10)
    suffix_length = width - prefix_length - 3
    prefix = ""
    suffix = ""

    for i in 1:length(str)
        if textwidth(str[1:prevind(str, i)]) >= prefix_length
            prefix = str[1:prevind(str, i)]
            break
        end
    end

    for i in length(str):-1:1
        if textwidth(str[nextind(str, i):end]) >= suffix_length
            suffix = str[nextind(str, i):end]
            break
        end
    end

    return "$(prefix)...$(suffix)"
end

_get_decomposition_data!(model::JuMP.Model) = get!(model.ext, :decomposition_jl, Dict{Symbol, Dict}())

function _print_iteration(k, args...)
    f(x) = Printf.@sprintf("%11.3e", x)
    println(lpad(k, 8), " │ ", join(f.(args), " │ "))
    return
end

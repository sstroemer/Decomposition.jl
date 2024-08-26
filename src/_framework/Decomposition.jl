module Decomposition

include("util/util.jl")

include("types/types.jl")

include("presolve/presolve.jl")
include("conversion/conversion.jl")
include("dualization/dualization.jl")
include("decomposition/decomposition.jl")

include("solve/solve.jl")

include("workflow/workflow.jl")

# This is directly taken from JuMP.jl and exports all internal symbols that do not start with an underscore (roughly).
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]
for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_") || startswith(sym_string, "@_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end])))
        continue
    end
    @eval export $sym
end

end

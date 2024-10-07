include("graphs.jl")
include("jump.jl")

function _print_iteration(k, args...)
    f(x) = (x isa Int) || (x isa AbstractString) ? lpad(x, 6) : Printf.@sprintf("%11.3e", x)
    output = string("│ ", f(k), " │ ", join(f.(args), " │ "), " │")
    println(output)
    return output
end

function _showtostr(obj::Any)
    sio = IOBuffer()
    show(sio, obj)
    return String(take!(sio))    
end

function Base.show(io::IO, attribute::AbstractDecompositionAttribute)
    str = "$(typeof(attribute))("
    str *= join(["$(prop)=$(getfield(attribute, prop))" for prop in propertynames(attribute)], ",")
    str *= ")"
    print(io, str)
    return nothing
end

function current_iteration(model::AbstractDecomposedModel)
    return length(model.info[:history])
end

gap_abs(x::Float64, y::Float64) = abs(x - y)

function gap_rel(x::Float64, y::Float64)
    tol = sqrt(eps(Float64))
    lower = min(x, y)
    upper = max(x, y)
    gap = upper - lower

    # This is similar to: https://www.gurobi.com/documentation/current/refman/mipgap2.html
    # Note: Behaviour for `upper == 0` may not be as expected.

    (isnan(gap) || !isfinite(gap)) && return +Inf
    isapprox(upper, 0.0; atol=tol) && return isapprox(lower, 0.0; atol=tol) ? 0.0 : +Inf
    return abs(gap / upper)
end

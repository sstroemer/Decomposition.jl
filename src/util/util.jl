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

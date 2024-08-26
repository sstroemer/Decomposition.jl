include("jump.jl")

dualize(node::String) = dualize(from_file(node))

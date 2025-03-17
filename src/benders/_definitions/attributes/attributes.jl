const _ATTRIBUTE_DEFAULTS = Set{AbstractGeneralAttribute}()

include("main.jl")
include("sub.jl")
include("config.jl")
include("termination.jl")

import .Main
import .Sub
import .Config
import .Termination

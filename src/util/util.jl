import SparseArrays
import UUIDs

import JuMP
const MOI = JuMP.MOI
import PaPILO


macro critical(msg, args...)
    msg = string(msg)
    return esc(quote
        @error $msg $(args...)
        error($msg)
    end)
end

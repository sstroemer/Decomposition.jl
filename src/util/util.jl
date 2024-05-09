import SparseArrays
import UUIDs
import Suppressor

import JuMP
import Dualization
const MOI = JuMP.MOI

import PaPILO
import SCIP_PaPILO_jll


macro critical(msg, args...)
    msg = string(msg)
    return esc(quote
        @error $msg $(args...)
        error($msg)
    end)
end

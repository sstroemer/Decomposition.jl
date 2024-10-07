@kwdef struct DualizeModel <: AbstractAttribute
    model::Any
    activate::Bool = true
end

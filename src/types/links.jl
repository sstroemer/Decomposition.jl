struct Link{T <: AbstractLink}
    Link{T}(from::AbstractNode, to::AbstractNode) where {T <: AbstractLink} = Link{T}(from, to, Dict{Symbol, Any}())
    Link{T}(from::AbstractNode, to::AbstractNode, content::Pair{Symbol, Any}) where {T <: AbstractLink} = Link{T}(from, to, Dict(content))
    function Link{T}(from::AbstractNode, to::AbstractNode, content::Dict{Symbol, Any}) where {T <: AbstractLink}
        link = T(from=from, to=to, info=content)
        push!(from._links, link)
        return link
    end
end

@kwdef struct DualLink <: AbstractLink
    from::AbstractNode
    to::AbstractNode
    info::Dict{Symbol, Any}
end

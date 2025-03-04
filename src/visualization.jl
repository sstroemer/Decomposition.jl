using PlotlyJS

function compare(hashes::Vector{String}, f_trace::Union{Function, Vector{Function}}, pltargs...; html::String, additional_traces = nothing, pltkwargs...)
    if f_trace isa Function
        f_trace = [f_trace]
    end
    
    path = normpath("out")

    # Iterate over hashes, instead of checking each file against the full list.
    # This allows controlling the order of traces, via the order of hashes.
    # A possible regex would instead be:
    #       `regex = Regex("_($(join(prefixes, "|")))\\.djl\\.json\$")`

    files = String[]
    for hash in hashes
        for elem in readdir(path)
            isfile(normpath(path, elem)) || continue
            endswith(elem, "_$(hash).djl.json") || continue
            push!(files, normpath(path, elem))
            break
        end
    end

    length(files) != length(hashes) && error("Could not find all files")
    
    traces = AbstractTrace[]
    !isnothing(additional_traces) && append!(traces, additional_traces)

    for filename in files
        model_info = JSON3.read(filename, Dict; allow_inf=true)
        push!(traces, scatter(; merge([f(model_info) for f in f_trace]...)...))
    end

    p = plot(traces, pltargs...; pltkwargs...)
    open("out/figures/$(html).html", "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
    return p
end

function trace_rel_gap(model_info::Dict)
    bounds = calculate_bounds(model_info)
    return (
        y=rel_gap.(bounds...),
        mode="lines",
        line_shape="hv",
    )
end

function trace_bounds(model_info::Dict)
    bounds = calculate_bounds(model_info)
    return (
        (y=bounds[1], mode="lines", line_shape="hv"),
        (y=bounds[2], mode="lines", line_shape="hv"),
    )
end

function calculate_bounds(model_info::Dict)
    lb = Vector{Float64}([model_info["history"][1]["lower_bound"]])
    ub = Vector{Float64}([model_info["history"][1]["upper_bound"]])
    for entry in model_info["history"][2:end]
        push!(lb, max(lb[end], entry["lower_bound"]))
        push!(ub, min(ub[end], entry["upper_bound"]))
    end
    return lb, ub
end

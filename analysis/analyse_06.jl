import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = replace(RUN_DIR, "experiments" => "analysis")
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

# iterations, sub time distribution
τ = 0.5
π0 = 1e8

y = Dict("h" => Dict{String, Any}("sub" => nothing), "g" => Dict{String, Any}("sub" => nothing))

# Extract results.
for r in RUNS
    for solver in ["h", "g"]
        dir = joinpath(RUN_DIR, r)
        timings = [(π = π0 * τ^(i - 1), t = JSON3.read(fn)) for (i, fn) in enumerate(readdir(dir; join=true)) if startswith(basename(fn), "$(solver)_timer")]

        y[solver]["x"] = [el.π for el in timings]
        y[solver]["iter"] = [el.t[:inner_timers][:main][:inner_timers][:optimize][:n_calls] for el in timings]

        tmp = []
        for timing in timings
            st = [v[:time_ns] for (k, v) in timing.t[:inner_timers][:sub][:inner_timers]]
            threshold = maximum(st) * 0.5
            st = mean(el for el in st if el > threshold) / 1e9
            push!(tmp, st)
        end
        y[solver]["sub"] = hcomb(y[solver]["sub"], tmp)
    end
end

base_iter = y["g"]["iter"][end] ./ 100.0
base_sub = y["g"]["sub"][end] ./ 100.0
for solver in ["h", "g"]
    # Average results.
    if length(RUNS) > 1
        y[solver]["sub"] = vec(mean(y[solver]["sub"]; dims = 2))
        y[solver]["sub"] = vec(mean(y[solver]["sub"]; dims = 2))
    end

    # Normalize.
    y[solver]["iter"] = y[solver]["iter"] ./ base_iter
    y[solver]["sub"] = y[solver]["sub"] ./ base_sub
end

# Plot.
function make_plot(traces, layout)
    kwlay = Dict(
        :title => "",
        :xaxis => PlotlyJS.attr(;
            showgrid = true,
            zeroline = false,
            showline = true,
            mirror = true,
            ticks = "outside",
            ticklen = 5,
            tickwidth = 1.5,
            tickcolor = "black",
            gridcolor = "lightgray",
        ),
        :yaxis => PlotlyJS.attr(;
            showgrid = true,
            zeroline = false,
            showline = true,
            mirror = true,
            ticks = "outside",
            ticklen = 5,
            tickwidth = 1.5,
            tickcolor = "black",
            gridcolor = "lightgray",
        ),
        :legend => PlotlyJS.attr(;
            x = 0.95,
            y = 0.95,
            bordercolor = "black",
            borderwidth = 1,
            xanchor = "right",
            yanchor = "top",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 12, color = "black"),
        :margin => PlotlyJS.attr(; l = 80, r = 50, b = 65, t = 90),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

x = reverse(string.(log.(1/τ, y["g"]["x"] ./ minimum(y["g"]["x"]))))

# Plot iterations.
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = reverse(y["g"]["iter"]), name = "Gurobi v12.0.1", marker_color = "#ff4800"))
push!(traces, bar(; x = x, y = reverse(y["h"]["iter"]), name = "HiGHS v1.9.0", marker_color = "#0026ff"))
savefig(
    make_plot(traces, (xaxis_title = "feasibility penalty (π ⋅ τ^x)", yaxis_title = "iterations (%)")),
    joinpath(VIZ_DIR, "iterations.png"), width = 500, height = 500,
)

# Plot avg. sub time.
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = reverse(y["g"]["sub"]), name = "Gurobi v12.0.1", marker_color = "#ff4800"))
push!(traces, bar(; x = x, y = reverse(y["h"]["sub"]), name = "HiGHS v1.9.0", marker_color = "#0026ff"))
savefig(
    make_plot(traces, (xaxis_title = "feasibility penalty (π ⋅ τ^x)", yaxis_title = "avg. sub-model solve time (%)")),
    joinpath(VIZ_DIR, "sub_time.png"), width = 500, height = 500,
)

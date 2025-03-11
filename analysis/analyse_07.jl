import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, bar, Layout, savefig

const EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
const RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
const RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
const RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

const examples = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "baseline"]
y = Dict(e => Dict{String, Any}("iter" => 0, "main" => 0.0, "main_aux" => 0.0, "sub" => 0.0) for e in examples)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for e in examples
        timings = JSON3.read(joinpath(dir, "timer_$(e).json"))
        mit = timings[:inner_timers][:main][:inner_timers]
        y[e]["iter"] += mit[:optimize][:n_calls]
        y[e]["main"] += mit[:optimize][:time_ns]
        y[e]["main_aux"] += sum(v[:time_ns] for (k, v) in mit if k != :optimize)
    end
end

# Average results (over all runs [already included in baseline], then down to "per iteration"), normalize to baseline.
baseline_iter = y["baseline"]["iter"] / 100.
baseline = (y["baseline"]["main"] + y["baseline"]["main_aux"]) / y["baseline"]["iter"] / 100.
for e in examples
    y[e]["main"] /= baseline * y[e]["iter"]
    y[e]["main_aux"] /= baseline * y[e]["iter"]
    y[e]["iter"] /= baseline_iter
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

exs = sort(examples; rev=true)
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = [y[e]["main_aux"] + y[e]["main"] for e in exs], y = exs, marker_color = "#0f4800", orientation="h", name="main-model (time overhead)", offsetgroup=1))
push!(traces, bar(; x = [y[e]["main"] for e in exs], y = exs, marker_color = "#ff4800", orientation="h", name="main-model (time solve)", offsetgroup=1))
push!(traces, bar(; x = [y[e]["iter"] for e in exs], y = exs, marker_color = "#0f48aa", orientation="h", name="main-model (iterations)", offsetgroup=2))
savefig(make_plot(traces, (barmode = "group", xaxis_title = "iterations / time compared to baseline (%)",)), joinpath(RUN_DIR, "fig.svg"), width = 900, height = 450)

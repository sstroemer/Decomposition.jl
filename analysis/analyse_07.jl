import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, bar, Layout, savefig

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = replace(RUN_DIR, "experiments" => "analysis")
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

examples = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "baseline"]
y = Dict(
    e => Dict{String, Any}("iter" => 0, "main" => 0.0, "main_aux" => 0.0, "sub" => 0.0, "sub_aux" => 0.0) for
    e in examples
)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for e in examples
        timings = JSON3.read(joinpath(dir, "timer_$(e).json"))

        mit = timings[:inner_timers][:main][:inner_timers]
        y[e]["iter"] += mit[:optimize][:n_calls]
        y[e]["main"] += mit[:optimize][:time_ns]
        y[e]["main_aux"] += sum(v[:time_ns] for (k, v) in mit if k != :optimize)

        sit = timings[:inner_timers][:sub][:inner_timers]
        worst = argmax(Dict(k => v[:time_ns] for (k, v) in sit))
        y[e]["sub"] += sit[worst][:inner_timers][:optimize][:time_ns]
        y[e]["sub_aux"] += sum(v[:time_ns] for (k, v) in sit[worst][:inner_timers] if k != :optimize)
    end
end

# Average results (over all runs [already included in baseline], then down to "per iteration"), normalize to baseline.
baseline_iter = y["baseline"]["iter"] / 100.0
baseline = sum(v for (k, v) in y["baseline"] if k != "iter") / y["baseline"]["iter"] / 100.0
for e in examples
    y[e]["main"] /= baseline * y[e]["iter"]
    y[e]["main_aux"] /= baseline * y[e]["iter"]
    y[e]["sub"] /= baseline * y[e]["iter"]
    y[e]["sub_aux"] /= baseline * y[e]["iter"]
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
            orientation = "h",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 12, color = "black"),
        :margin => PlotlyJS.attr(; l = 80, r = 50, b = 65, t = 90),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

names = Dict(
    "baseline" => "opt. cuts only",
    "ex1" => "opt. & feas. cuts",
    "ex2" => "opt. & feas. + merge small",
    "ex3" => "opt. & feas. + merge all",
    "ex4" => "opt. & feas. cuts + 1sub",
    "ex5" => "opt. & feas. + merge small + 1sub",
    "ex6" => "opt. & feas. + merge all + 1sub",
)

exs = sort(examples; rev = true)
yn = [names[e] for e in exs]

traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    bar(;
        x = [y[e]["main_aux"] + y[e]["main"] + y[e]["sub_aux"] + y[e]["sub"] for e in exs],
        y = yn,
        marker_color = "#7ea15c",
        orientation = "h",
        name = "time (overhead)",
        offsetgroup = 1,
        legendgrouptitle = PlotlyJS.attr(; text = "sub (worst)"),
        legendgroup = "sub",
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main_aux"] + y[e]["main"] + y[e]["sub"] for e in exs],
        y = yn,
        marker_color = "#458a00",
        orientation = "h",
        name = "time (solve)",
        legendgrouptitle = PlotlyJS.attr(; text = "sub (worst)"),
        legendgroup = "sub",
        offsetgroup = 1,
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main_aux"] + y[e]["main"] for e in exs],
        y = yn,
        marker_color = "#b85c5c",
        orientation = "h",
        name = "time (overhead)",
        offsetgroup = 1,
        legendgrouptitle = PlotlyJS.attr(; text = "main"),
        legendgroup = "main",
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main"] for e in exs],
        y = yn,
        marker_color = "#b80000",
        orientation = "h",
        name = "time (solve)",
        legendgrouptitle = PlotlyJS.attr(; text = "main"),
        legendgroup = "main",
        offsetgroup = 1,
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["iter"] for e in exs],
        y = yn,
        marker_color = "#0f48aa",
        orientation = "h",
        name = "iterations",
        legendgrouptitle = PlotlyJS.attr(; text = "general"),
        legendgroup = "general",
        offsetgroup = 2,
    ),
)
savefig(
    make_plot(traces, (barmode = "group", xaxis_title = "iterations / time compared to baseline (%)")),
    joinpath(VIZ_DIR, "fig.png");
    width = 900,
    height = 400,
)

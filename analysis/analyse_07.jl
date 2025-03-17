import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, bar, Layout, savefig

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

examples = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "baseline"]
y = Dict(e => Dict{String, Any}("iter" => 0, "main" => 0.0, "sub_p" => 0.0, "sub_s" => 0.0) for e in examples)

PARALLELIZATION = 16

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for e in examples
        timings = JSON3.read(joinpath(dir, "timer_$(e).json"))
        mit = timings[:inner_timers][:main]
        if mit[:inner_timers][:optimize][:n_calls] >= 250
            @error "not converged" r dir e
            continue
        end

        y[e]["iter"] += mit[:inner_timers][:optimize][:n_calls]
        y[e]["main"] += mit[:time_ns]

        sit = timings[:inner_timers][:sub][:inner_timers]
        par = zeros(PARALLELIZATION)
        pi = 1
        for t in sort([v[:time_ns] for v in values(sit)]; rev = true)
            par[pi] += t
            pi = pi % PARALLELIZATION + 1
        end
        y[e]["sub_p"] += par[1]
        y[e]["sub_s"] += sum(v[:time_ns] for (k, v) in sit)
    end
end

# Average results (over all runs [already included in baseline], then down to "per iteration"), normalize to baseline.
baseline_iter = y["baseline"]["iter"] / 100.0
baseline = (y["baseline"]["main"] + y["baseline"]["sub_p"]) / 100.0
for e in examples
    y[e]["main"] /= baseline
    y[e]["sub_p"] /= baseline
    y[e]["sub_s"] /= baseline
    y[e]["iter"] /= baseline_iter
end

# Plot.
function make_plot(traces, layout)
    kwlay = Dict(
        :title => "",
        :titlefont_size => 12,
        :xaxis_type => "log",
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
            x = 1.0,
            y = 0.5,
            bordercolor = "black",
            borderwidth = 1,
            xanchor = "left",
            yanchor = "center",
            orientation = "v",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 12, color = "black"),
        :margin => PlotlyJS.attr(; l = 60, r = 50, b = 65, t = 40),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

names = Dict(
    "baseline" => "optimality cuts only",
    "ex1" => "opt. & feas. cuts",
    "ex2" => "o./f. cuts & small-sub VIs",
    "ex4" => "o./f. cuts & 1st-sub VIs",
    "ex5" => "o./f. cuts & small-and-1st-sub VIs",
)

exs = sort(collect(e for e in examples if haskey(names, e)); rev = true)
yn = [names[e] for e in exs]

traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    bar(;
        x = [y[e]["iter"] for e in exs],
        y = yn,
        marker_color = "#0f48aa",
        orientation = "h",
        name = "iterations",
        offsetgroup = 1,
        # legendgrouptitle = PlotlyJS.attr(; text = "time"),
        # legendgroup = "iterations",
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main"] + y[e]["sub_s"] for e in exs],
        y = yn,
        marker_color = "#7ea15c",
        orientation = "h",
        name = "sub (serial)",
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "model time"),
        legendgroup = "time",
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main"] + y[e]["sub_p"] for e in exs],
        y = yn,
        marker_color = "#458a00",
        orientation = "h",
        name = "sub (parallel)",
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "model time"),
        legendgroup = "time",
    ),
)
push!(
    traces,
    bar(;
        x = [y[e]["main"] for e in exs],
        y = yn,
        marker_color = "#b85c5c",
        orientation = "h",
        name = "main",
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "time"),
        legendgroup = "time",
    ),
)
savefig(
    make_plot(traces, (barmode = "group", xaxis_title = "iterations / time compared to baseline (%)")),
    joinpath(VIZ_DIR, "fig.png");
    width = 900,
    height = 400,
)

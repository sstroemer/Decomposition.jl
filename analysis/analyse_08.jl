import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, scatter, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

exp_tol = collect(0:3)
tol = vcat([0.0], 10.0 .^ (-exp_tol[2:end]))

y = Dict(i => Dict{String, Any}("iter_1" => 0, "time_1" => 0, "iter_5" => 0, "time_5" => 0) for i in exp_tol)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for i in exp_tol
        history = JSON3.read(joinpath(dir, "history_$(i).json"))
        # timings = JSON3.read(joinpath(dir, "timer_$(i).json"))

        t0 = history[1]["t"]["timestamp"]

        for entry in history
            gap = (entry["ub"] - entry["lb"]) / entry["ub"]
            if gap <= 0.01
                y[i]["iter_1"] += entry["k"]
                y[i]["time_1"] += entry["t"]["timestamp"] - t0
                break
            end
        end

        for entry in history
            gap = (entry["ub"] - entry["lb"]) / entry["ub"]
            if gap <= 0.05
                y[i]["iter_5"] += entry["k"]
                y[i]["time_5"] += entry["t"]["timestamp"] - t0
                break
            end
        end
    end
end

# Normalize.
base_iter = y[0]["iter_1"] / 100.0
base_time = y[0]["time_1"] / 100.0
for i in exp_tol
    y[i]["iter_1"] /= base_iter
    y[i]["time_1"] /= base_time
    y[i]["iter_5"] /= base_iter
    y[i]["time_5"] /= base_time
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
            x = 0.50,
            y = 0.95,
            bordercolor = "black",
            borderwidth = 1,
            xanchor = "center",
            yanchor = "bottom",
            orientation = "h",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 10, color = "black"),
        :margin => PlotlyJS.attr(; l = 0, r = 10, b = 10, t = 10),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

x = vcat(["default (1e-6)"], ["1e-$(i)" for i in exp_tol[2:end]])
ymax = maximum(maximum.(values.(values(y))))

# Plot iterations.
traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    bar(;
        x = x,
        y = vcat(y[0]["iter_1"], zeros(length(exp_tol) - 1)),
        offsetgroup = 1,
        legendgroup = "concurrent",
        name = "1% gap",
        marker_color = "#9c0000ff",
        legendgrouptitle = PlotlyJS.attr(; text = "concurrent with crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(y[0]["iter_5"], zeros(length(exp_tol) - 1)),
        offsetgroup = 2,
        legendgroup = "concurrent",
        name = "5% gap",
        marker_color = "#9c000077",
        legendgrouptitle = PlotlyJS.attr(; text = "concurrent with crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(0.0, [y[i]["iter_1"] for i in exp_tol[2:end]]),
        offsetgroup = 1,
        legendgroup = "barrier",
        name = "1% gap",
        marker_color = "#00009cff",
        legendgrouptitle = PlotlyJS.attr(; text = "barrier w/o crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(0.0, [y[i]["iter_5"] for i in exp_tol[2:end]]),
        offsetgroup = 2,
        legendgroup = "barrier",
        name = "5% gap",
        marker_color = "#00009c77",
        legendgrouptitle = PlotlyJS.attr(; text = "barrier w/o crossover"),
    ),
)
savefig(
    make_plot(
        traces,
        (xaxis_title = "barrier convergence tolerance", yaxis_title = "iterations (%)", yaxis_range = [0, ymax * 1.1]),
    ),
    joinpath(VIZ_DIR, "iterations.png");
    width = 500,
    height = 500,
)

# Plot time.
traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    bar(;
        x = x,
        y = vcat(y[0]["time_1"], zeros(length(exp_tol) - 1)),
        offsetgroup = 1,
        legendgroup = "concurrent",
        name = "1% gap",
        marker_color = "#9c0000ff",
        legendgrouptitle = PlotlyJS.attr(; text = "concurrent with crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(y[0]["time_5"], zeros(length(exp_tol) - 1)),
        offsetgroup = 2,
        legendgroup = "concurrent",
        name = "5% gap",
        marker_color = "#9c000077",
        legendgrouptitle = PlotlyJS.attr(; text = "concurrent with crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(0.0, [y[i]["time_1"] for i in exp_tol[2:end]]),
        offsetgroup = 1,
        legendgroup = "barrier",
        name = "1% gap",
        marker_color = "#00009cff",
        legendgrouptitle = PlotlyJS.attr(; text = "barrier w/o crossover"),
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = vcat(0.0, [y[i]["time_5"] for i in exp_tol[2:end]]),
        offsetgroup = 2,
        legendgroup = "barrier",
        name = "5% gap",
        marker_color = "#00009c77",
        legendgrouptitle = PlotlyJS.attr(; text = "barrier w/o crossover"),
    ),
)
savefig(
    make_plot(
        traces,
        (xaxis_title = "barrier convergence tolerance", yaxis_title = "time (%)", yaxis_range = [0, ymax * 1.1]),
    ),
    joinpath(VIZ_DIR, "time.png");
    width = 500,
    height = 500,
)

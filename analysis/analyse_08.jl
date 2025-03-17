import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, scatter, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = replace(RUN_DIR, "experiments" => "analysis")
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

exp_tol = collect(0:10)
tol = vcat([0.0], 10.0 .^ (-exp_tol[2:end]))

y = Dict(
    i => Dict{String, Any}("iter" => nothing, "time" => nothing, "lb" => nothing, "ub" => nothing) for i in exp_tol
)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for i in exp_tol
        history = JSON3.read(joinpath(dir, "history_$(i).json"))
        # timings = JSON3.read(joinpath(dir, "timer_$(i).json"))

        y[i]["iter"] = [h["k"] for h in history]

        t0 = history[1]["t"]["timestamp"]
        y[i]["time"] = hcomb(y[i]["time"], [(h["t"]["timestamp"] - t0) for h in history])

        lb, ub = [-Inf], [+Inf]
        for j in eachindex(history)
            push!(lb, max(lb[end], history[j]["lb"]))
            push!(ub, min(ub[end], history[j]["ub"]))
        end
        y[i]["lb"] = hcomb(y[i]["lb"], lb[2:end])
        y[i]["ub"] = hcomb(y[i]["ub"], ub[2:end])
    end
end

# Average results and convert to seconds.
for i in exp_tol
    y[i]["time"] = vec(mean(y[i]["time"]; dims = 2)) ./ 1e9
    y[i]["lb"] = vec(mean(y[i]["lb"]; dims = 2))
    y[i]["ub"] = vec(mean(y[i]["ub"]; dims = 2))
end

# Precalculate relative gap.
for i in exp_tol
    y[i]["gap"] = (y[i]["ub"] - y[i]["lb"]) ./ y[i]["ub"]
    y[i]["gap_reached"] = [findfirst(<=(th), y[i]["gap"]) for th in [0.1, 0.05, 0.025, 0.01]]
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
            x = 0.10,
            y = 0.95,
            bordercolor = "black",
            borderwidth = 1,
            xanchor = "left",
            yanchor = "top",
            orientation = "h",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 10, color = "black"),
        :margin => PlotlyJS.attr(; l = 0, r = 10, b = 10, t = 10),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

opacity = ["ff", "cf", "b0", "9f"]

iter_max = maximum(y[i]["iter"][end] for i in exp_tol)
time_max = maximum(y[i]["time"][end] for i in exp_tol)
base_iter = y[0]["iter"][end] / 100.0
base_time = y[0]["time"][end] / 100.0
ymax = max(iter_max / base_iter, time_max / base_time) * 1.05
x = [(e == 0) ? "" : "1e-$(e)" for e in exp_tol]

traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    scatter(;
        zorder = -1,
        x = x,
        y = fill(100.0, length(x)),
        mode = "lines",
        showlegend = false,
        line = PlotlyJS.attr(; color = "black", dash = "dash"),
    ),
)
for (i, gr) in enumerate(y[0]["gap_reached"])
    push!(
        traces,
        bar(;
            zorder = -2,
            x = x,
            y = vcat(y[0]["iter"][gr] / base_iter, zeros(length(exp_tol) - 1)),
            marker_color = "#ff4800$(opacity[i])",
            legendgrouptitle = PlotlyJS.attr(; text = "concurrent"),
            legendgroup = "concurrent",
            name = "Δ = $(["10%", "5%", "2.5%", "1%"][i])",
        ),
    )
    push!(
        traces,
        bar(;
            zorder = 1,
            x = x,
            y = vcat(0.0, [y[e]["iter"][y[e]["gap_reached"][i]] / base_iter for e in exp_tol[2:end]]),
            marker_color = "#0026ff$(opacity[i])",
            legendgrouptitle = PlotlyJS.attr(; text = "barrier"),
            legendgroup = "barrier",
            name = "Δ = $(["10%", "5%", "2.5%", "1%"][i])",
        ),
    )
end
savefig(
    make_plot(
        traces,
        (
            barmode = "overlay",
            xaxis_title = "barrier tolerance",
            yaxis_title = "iterations (%)",
            yaxis_range = [0, ymax],
        ),
    ),
    joinpath(VIZ_DIR, "iterations.png");
    width = 700,
    height = 550,
)

traces = Vector{PlotlyJS.GenericTrace}()
push!(
    traces,
    scatter(;
        zorder = -1,
        x = x,
        y = fill(100.0, length(x)),
        mode = "lines",
        showlegend = false,
        line = PlotlyJS.attr(; color = "black", dash = "dash"),
    ),
)
for (i, gr) in enumerate(y[0]["gap_reached"])
    push!(
        traces,
        bar(;
            zorder = -2,
            x = x,
            y = vcat(y[0]["time"][gr] / base_time, zeros(length(exp_tol) - 1)),
            marker_color = "#ff4800$(opacity[i])",
            legendgrouptitle = PlotlyJS.attr(; text = "concurrent"),
            legendgroup = "concurrent",
            name = "Δ = $(["10%", "5%", "2.5%", "1%"][i])",
        ),
    )
    push!(
        traces,
        bar(;
            zorder = 1,
            x = x,
            y = vcat(0.0, [y[e]["time"][y[e]["gap_reached"][i]] / base_time for e in exp_tol[2:end]]),
            marker_color = "#0026ff$(opacity[i])",
            legendgrouptitle = PlotlyJS.attr(; text = "barrier"),
            legendgroup = "barrier",
            name = "Δ = $(["10%", "5%", "2.5%", "1%"][i])",
        ),
    )
end
savefig(
    make_plot(
        traces,
        (barmode = "overlay", xaxis_title = "barrier tolerance", yaxis_title = "time (%)", yaxis_range = [0, ymax]),
    ),
    joinpath(VIZ_DIR, "time.png");
    width = 700,
    height = 550,
)

savefig(
    make_plot(
        traces,
        (barmode = "overlay", xaxis_title = "barrier tolerance", yaxis_title = "time (%)", yaxis_range = [0, ymax]),
    ),
    joinpath(VIZ_DIR, "time.png");
    scale = 2,
    width = 700,
    height = 550,
)

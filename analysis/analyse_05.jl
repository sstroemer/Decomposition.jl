import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, scatter, Layout, savefig

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = replace(RUN_DIR, "experiments" => "analysis")
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

examples = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6"]
y = Dict(e => Dict{String, Any}("lb" => nothing, "ub" => nothing) for e in examples)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for e in examples
        history = JSON3.read(joinpath(dir, "history_$(e).json"))

        lb, ub = [-Inf], [+Inf]
        for i in eachindex(history)
            push!(lb, max(lb[end], history[i]["lb"]))
            push!(ub, min(ub[end], history[i]["ub"]))
        end

        y[e]["lb"] = hcomb(y[e]["lb"], lb[2:end])
        y[e]["ub"] = hcomb(y[e]["ub"], ub[2:end])
    end
end

# Average results.
if length(RUNS) > 1
    for e in examples
        y[e]["lb"] = vec(mean(y[e]["lb"]; dims = 2))
        y[e]["ub"] = vec(mean(y[e]["ub"]; dims = 2))
    end
end

# Get min/max bounds.
ymin = minimum(minimum(y[e]["lb"]) for e in examples)
ymax = maximum(maximum(y[e]["ub"]) for e in examples)
xmin = 1
xmax = maximum(length(y[e]["lb"]) for e in examples)

# Fill with average after termination due to convergence.
for e in examples
    if length(y[e]["lb"]) < xmax
        avg = (y[e]["lb"][end] + y[e]["ub"][end]) / 2
        filler = fill(avg, xmax - length(y[e]["lb"]))
        y[e]["lb"] = vcat(y[e]["lb"], filler)
        y[e]["ub"] = vcat(y[e]["ub"], filler)
    end
end

# Plot.
function make_plot(ex)
    x = xmin:xmax
    colors = ["#0026ff", "#ff4800", "#ffa480"]
    legend_entries = Dict(
        "ex1" => "baseline",
        "ex2" => "bounded variables",
        "ex3" => "baseline",
        "ex4" => "bounded variables",
        "ex5" => "baseline",
        "ex6" => "bounded variables",
    )

    traces = Vector{PlotlyJS.GenericTrace}()
    push!(
        traces,
        scatter(;
            x = x,
            y = fill(5.821398080e+06, length(x)),
            mode = "lines",
            name = "target",
            line = PlotlyJS.attr(; color = "black", dash = "dash"),
        ),
    )
    for (i, e) in enumerate(ex)
        color = colors[i]
        x = xmin:xmax
        push!(traces, scatter(; x = x, y = y[e]["lb"], mode = "lines", name = legend_entries[e], line_color = color))
        push!(traces, scatter(; x = x, y = y[e]["ub"], mode = "lines", showlegend = false, line_color = color))
    end

    return plot(
        traces,
        Layout(;
            title = "",
            xaxis_title = "iteration",
            yaxis_title = "best objective bounds",
            yaxis_type = "log",
            xaxis = PlotlyJS.attr(;
                range = [1, xmax],
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
            yaxis = PlotlyJS.attr(;
                range = [0, log10(ymax) * 1.05],
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
            legend = PlotlyJS.attr(;
                x = 0.95,
                y = 0.95,
                bordercolor = "black",
                borderwidth = 1,
                xanchor = "right",
                yanchor = "top",
            ),
            plot_bgcolor = "white",
            paper_bgcolor = "white",
            font = PlotlyJS.attr(; family = "Arial, sans-serif", size = 12, color = "black"),
            margin = PlotlyJS.attr(; l = 80, r = 50, b = 65, t = 90),
        ),
    )
end

savefig(make_plot(examples[1:2]), joinpath(VIZ_DIR, "simplex.png"); width = 450, height = 550)
savefig(make_plot(examples[3:4]), joinpath(VIZ_DIR, "ipm.png"); width = 450, height = 550)
savefig(make_plot(examples[5:6]), joinpath(VIZ_DIR, "stabilized.png"); width = 450, height = 550)

import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, scatter, Layout, savefig

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

COLORS = ["#9c0000ff", "#009c00ff", "#00009cff"]
COLORS_T = ["#9c000077", "#009c0077", "#00009c77"]

y = Dict()
y_keys = ["lb", "ub", "t"]
xmax = 250

files = []
for r in RUNS
    fh, ft = readdir(joinpath(RUN_DIR, r); join = true)
    i = parse(Int64, rsplit(rsplit(fh, "."; limit = 2)[1], "_"; limit = 2)[2])
    push!(files, (i, fh, ft))
end

for f in files
    n = f[1]
    haskey(y, n) || (y[n] = Dict{String, Any}(k => nothing for k in y_keys))
    history = JSON3.read(f[2]; allow_inf = true)
    timings = JSON3.read(f[3]; allow_inf = true)

    t0 = history[1]["t"]["timestamp"]

    lb, ub, t = [-Inf], [+Inf], [0.0]
    for i in eachindex(history)
        push!(lb, max(lb[end], history[i]["lb"]))
        push!(ub, min(ub[end], abs(history[i]["ub"])))
        push!(t, history[i]["t"]["timestamp"] - t0)
    end
    for i in 1:(xmax-length(history))
        push!(lb, lb[end])
        push!(ub, ub[end])
        push!(t, t[end])
    end

    y[n]["lb"] = hcomb(y[n]["lb"], lb[2:(xmax+1)])
    y[n]["ub"] = hcomb(y[n]["ub"], ub[2:(xmax+1)])
    y[n]["t"] = hcomb(y[n]["t"], t[2:(xmax+1)])
end

@info "EXPERIMENT $(EXPERIMENT_NR)" nof_experiments = length(y) avg_runs = length(files) / length(y)

# Average results.
for e in keys(y)
    y[e]["lb"] isa Vector && continue
    y[e]["lb"] = vec(mean(y[e]["lb"]; dims = 2))
    y[e]["ub"] = vec(mean(y[e]["ub"]; dims = 2))
    y[e]["t"] = vec(mean(y[e]["t"]; dims = 2))
end

# Prepare per iteration timings.
tmax = maximum([e["t"][end] for e in values(y)])
tbase = y[1]["t"][end] / tmax
iterbase = findfirst(el -> el <= 1e-2, (y[1]["ub"] .- y[1]["lb"]) ./ y[1]["ub"])
yt = Dict()
for e in keys(y)
    xt = 1:1000
    yt[e] = Dict("t" => zeros(1000), "lb" => zeros(1000), "ub" => zeros(1000))
    for t in xt
        yt[e]["t"][t] = t / 1000.0 * 100.0
        idx = findlast(x -> x <= t / 1000.0 * tmax, y[e]["t"])
        yt[e]["lb"][t] = y[e]["lb"][idx]
        yt[e]["ub"][t] = y[e]["ub"][idx]
    end
end

# Fill with average after termination due to convergence.
# for e in examples
#     if length(y[e]["lb"]) < xmax
#         avg = (y[e]["lb"][end] + y[e]["ub"][end]) / 2
#         filler = fill(avg, xmax - length(y[e]["lb"]))
#         y[e]["lb"] = vcat(y[e]["lb"], filler)
#         y[e]["ub"] = vcat(y[e]["ub"], filler)
#     end
# end

# Plot.
function make_plot(cx, cy, ex, xaxt)
    colors = [COLORS[1], COLORS[3]]
    legend_entries = Dict(
        1 => "baseline",
        2 => "bounded variables",
        3 => "baseline",
        4 => "bounded variables",
        5 => "baseline",
        6 => "bounded variables",
    )

    traces = Vector{PlotlyJS.GenericTrace}()
    annotations = []
    ax = 0

    push!(
        traces,
        scatter(;
            x = cx,
            y = fill(2.23922483671398e7, length(cx)),
            mode = "lines",
            name = "target",
            line = PlotlyJS.attr(; color = "black", dash = "dash"),
        ),
    )
    for (i, e) in enumerate(ex)
        color = colors[i]
        push!(traces, scatter(; x = cx, y = cy[e]["lb"], mode = "lines", name = legend_entries[e], line_color = color))
        push!(traces, scatter(; x = cx, y = cy[e]["ub"], mode = "lines", showlegend = false, line_color = color))

        gap = (cy[e]["ub"] .- cy[e]["lb"]) ./ cy[e]["ub"]
        cidx = findfirst(gap .<= 1e-2)

        isnothing(cidx) && continue

        cidxy = log10((cy[e]["ub"][cidx] .+ cy[e]["lb"][cidx]) / 2)

        push!(
            annotations,
            PlotlyJS.attr(;
                x = cx[cidx],
                y = cidxy,
                ax = 0,
                ay = 0,
                text = "<b>×</b>",
                arrowcolor = color,
                font = PlotlyJS.attr(; color = color, size = 20),
            ),
        )
        ax += 1
    end

    return plot(
        traces,
        Layout(;
            annotations = annotations,
            title = "",
            xaxis_title = "% of maximum $(xaxt)",
            yaxis_title = "best objective bounds",
            yaxis_type = "log",
            xaxis = PlotlyJS.attr(;
                range = [1, cx[end]],
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
                range = [2, 14],
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

savefig(
    make_plot((1:250) ./ iterbase .* 100, y, 1:2, "iterations"),
    joinpath(VIZ_DIR, "iter_0_simplex.png");
    width = 450,
    height = 550,
)
savefig(
    make_plot((1:250) ./ iterbase .* 100, y, 3:4, "iterations"),
    joinpath(VIZ_DIR, "iter_1_ipm.png");
    width = 450,
    height = 550,
)
savefig(
    make_plot((1:250) ./ iterbase .* 100, y, 5:6, "iterations"),
    joinpath(VIZ_DIR, "iter_2_stabilized.png");
    width = 450,
    height = 550,
)

savefig(
    make_plot((1:1000) ./ tbase ./ 10.0, yt, 1:2, "time"),
    joinpath(VIZ_DIR, "time_0_simplex.png");
    width = 450,
    height = 550,
)
savefig(
    make_plot((1:1000) ./ tbase ./ 10.0, yt, 3:4, "time"),
    joinpath(VIZ_DIR, "time_1_ipm.png");
    width = 450,
    height = 550,
)
savefig(
    make_plot((1:1000) ./ tbase ./ 10.0, yt, 5:6, "time"),
    joinpath(VIZ_DIR, "time_2_stabilized.png");
    width = 450,
    height = 550,
)

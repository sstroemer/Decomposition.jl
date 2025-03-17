import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

# iterations, sub time distribution
τ = 0.5
π0 = 1e8

PARALLELIZATION = 16
y = Dict(s => Dict{String, Any}(k => nothing for k in ["x", "iter", "main", "sub_p", "sub_s"]) for s in ["h", "g"])

# Extract results.
for r in RUNS
    for solver in ["h", "g"]
        dir = joinpath(RUN_DIR, r)
        timings = [
            (π = π0 * τ^(i - 1), t = JSON3.read(fn)) for (i, fn) in
            enumerate([fn for fn in readdir(dir; join = true) if startswith(basename(fn), "$(solver)_timer")])
        ]

        y[solver]["x"] = hcomb(y[solver]["x"], [el.π for el in timings])
        y[solver]["iter"] = hcomb(y[solver]["iter"], [el.t[:inner_timers][:sub][:n_calls] for el in timings])
        y[solver]["main"] = hcomb(y[solver]["main"], [el.t[:inner_timers][:main][:time_ns] for el in timings])
        y[solver]["sub_s"] = hcomb(y[solver]["sub_s"], [el.t[:inner_timers][:sub][:time_ns] for el in timings])

        tmp = []
        for el in timings
            sit = el.t[:inner_timers][:sub][:inner_timers]
            par = zeros(PARALLELIZATION)
            pi = 1
            for t in sort([v[:time_ns] for v in values(sit)]; rev = true)
                par[pi] += t
                pi = pi % PARALLELIZATION + 1
            end
            push!(tmp, par[1])
        end
        y[solver]["sub_p"] = hcomb(y[solver]["sub_p"], tmp)
    end
end

# Average results.
for solver in ["h", "g"]
    if length(RUNS) > 1
        for k in keys(y[solver])
            y[solver][k] = vec(mean(y[solver][k]; dims = 2))
        end
    end
end

# Normalize.
base_iter = y["g"]["iter"][end] ./ 100.0
base_time = (y["g"]["main"][end] + y["g"]["sub_p"][end]) ./ 100.0
for solver in ["h", "g"]
    y[solver]["iter"] = y[solver]["iter"] ./ base_iter
    y[solver]["main"] = y[solver]["main"] ./ base_time
    y[solver]["sub_p"] = y[solver]["sub_p"] ./ base_time
    y[solver]["sub_s"] = y[solver]["sub_s"] ./ base_time
end

# Plot.
function make_plot(traces, layout)
    kwlay = Dict(
        :title => "",
        # :yaxis_type => "log",
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
            x = 1.1,
            y = 1.1,
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

x = reverse(string.(y["g"]["x"]))

# Plot iterations.
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = reverse(y["g"]["iter"]), name = "Gurobi v12.0.1", marker_color = "#9c0000"))
push!(traces, bar(; x = x, y = reverse(y["h"]["iter"]), name = "HiGHS v1.9.0", marker_color = "#39009c"))
savefig(
    make_plot(traces, (xaxis_title = "feasibility penalty", yaxis_title = "iterations (%)")),
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
        y = reverse(y["g"]["main"] .+ y["g"]["sub_s"]),
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "Gurobi v12.0.1"),
        legendgroup = "gurobi",
        name = "sub (serial)",
        marker_color = "#9c7a68",
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = reverse(y["g"]["main"] .+ y["g"]["sub_p"]),
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "Gurobi v12.0.1"),
        legendgroup = "gurobi",
        name = "sub (parallel)",
        marker_color = "#9c3600",
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = reverse(y["g"]["main"]),
        offsetgroup = 2,
        legendgrouptitle = PlotlyJS.attr(; text = "Gurobi v12.0.1"),
        legendgroup = "gurobi",
        name = "main",
        marker_color = "#9c0000",
    ),
)

push!(
    traces,
    bar(;
        x = x,
        y = reverse(y["h"]["main"] .+ y["h"]["sub_s"]),
        offsetgroup = 1,
        legendgrouptitle = PlotlyJS.attr(; text = "HiGHS v1.9.0"),
        legendgroup = "highs",
        name = "sub (serial)",
        marker_color = "#8d689c",
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = reverse(y["h"]["main"] .+ y["h"]["sub_p"]),
        offsetgroup = 1,
        legendgrouptitle = PlotlyJS.attr(; text = "HiGHS v1.9.0"),
        legendgroup = "highs",
        name = "sub (parallel)",
        marker_color = "#6f009c",
    ),
)
push!(
    traces,
    bar(;
        x = x,
        y = reverse(y["h"]["main"]),
        offsetgroup = 1,
        legendgrouptitle = PlotlyJS.attr(; text = "HiGHS v1.9.0"),
        legendgroup = "highs",
        name = "main",
        marker_color = "#39009c",
    ),
)

savefig(
    make_plot(traces, (xaxis_title = "feasibility penalty", yaxis_title = "total model time (%)")),
    joinpath(VIZ_DIR, "time.png");
    width = 500,
    height = 500,
)

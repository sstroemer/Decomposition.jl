import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))
hcomb(a, b) = isnothing(a) ? b : hcat(a, b)

PARALLELIZATION = 16
y = Dict()
y_keys = ["iter", "main", "sub_p", "sub_s", "cnt"]

files = []
for r in RUNS
    for ft in readdir(joinpath(RUN_DIR, r); join = true)
        tmp = rsplit(rsplit(ft, "."; limit = 2)[1], "_"; limit = 4)
        s, i = string(tmp[end-2][end]), parse(Int64, tmp[end])
        push!(files, (ft, s, i))
        y[s] = Dict()
    end
end

for f in files
    s, i = f[2:3]
    haskey(y[s], i) || (y[s][i] = Dict{String, Any}(k => 0 for k in y_keys))
    timings = JSON3.read(f[1]; allow_inf = true)

    y[s][i]["iter"] += timings[:inner_timers][:sub][:n_calls]
    y[s][i]["main"] += timings[:inner_timers][:main][:time_ns]
    y[s][i]["sub_s"] += timings[:inner_timers][:sub][:time_ns]

    sit = timings[:inner_timers][:sub][:inner_timers]
    par = zeros(PARALLELIZATION)
    pi = 1
    for t in sort([v[:time_ns] for v in values(sit)]; rev = true)
        par[pi] += t
        pi = pi % PARALLELIZATION + 1
    end
    y[s][i]["sub_p"] += par[1]

    y[s][i]["cnt"] += 1
end

@info "EXPERIMENT $(EXPERIMENT_NR)" nof_experiments_g = length(y["g"]) avg_runs_g =
    length(collect(f for f in files if contains(f[1], "g_timer"))) / length(y["g"])
@info "EXPERIMENT $(EXPERIMENT_NR)" nof_experiments_h = length(y["h"]) avg_runs_h =
    length(collect(f for f in files if contains(f[1], "h_timer"))) / length(y["h"])

# Average results.
for solver in ["h", "g"]
    for i in keys(y[solver])
        for k in keys(y[solver][i])
            (k == "cnt") && continue
            y[solver][i][k] /= y[solver][i]["cnt"]
        end
    end
end

# Normalize.
base_iter = y["g"][5]["iter"][end] ./ 100.0
base_time = (y["g"][5]["main"][end] + y["g"][5]["sub_p"][end]) ./ 100.0
for solver in ["h", "g"]
    for i in keys(y[solver])
        y[solver][i]["iter"] = y[solver][i]["iter"] ./ base_iter
        y[solver][i]["main"] = y[solver][i]["main"] ./ base_time
        y[solver][i]["sub_p"] = y[solver][i]["sub_p"] ./ base_time
        y[solver][i]["sub_s"] = y[solver][i]["sub_s"] ./ base_time
    end
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

x = sort(collect(keys(y["g"])))
xs = ["1e$(i)" for i in x]

# Plot iterations.
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = xs, y = [y["g"][i]["iter"] for i in x], name = "Gurobi v12.0.1", marker_color = "#9c0000"))
push!(traces, bar(; x = xs, y = [y["h"][i]["iter"] for i in x], name = "HiGHS v1.9.0", marker_color = "#39009c"))
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
        x = xs,
        y = [y["g"][i]["main"] + y["g"][i]["sub_s"] for i in x],
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
        x = xs,
        y = [y["g"][i]["main"] + y["g"][i]["sub_p"] for i in x],
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
        x = xs,
        y = [y["g"][i]["main"] for i in x],
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
        x = xs,
        y = [y["h"][i]["main"] + y["h"][i]["sub_s"] for i in x],
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
        x = xs,
        y = [y["h"][i]["main"] + y["h"][i]["sub_p"] for i in x],
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
        x = xs,
        y = [y["h"][i]["main"] for i in x],
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

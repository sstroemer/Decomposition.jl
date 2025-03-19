import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, bar, Layout, savefig

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
    ft = only(readdir(joinpath(RUN_DIR, r); join = true))
    n, drop, ipm = parse.([Int64, Int64, Bool], rsplit(rsplit(ft, "."; limit = 2)[1], "_"; limit = 4)[end-2:end])
    push!(files, (ft, n, drop, ipm))
end

for f in files
    e = n, drop, ipm = f[2:end]
    haskey(y, e) || (y[e] = Dict{String, Any}(k => 0.0 for k in y_keys))
    timings = JSON3.read(f[1]; allow_inf = true)

    mit = timings[:inner_timers][:main]
    if mit[:inner_timers][:optimize][:n_calls] >= 500
        @error "not converged" n drop ipm
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

    y[e]["cnt"] += 1
end

examples = collect(keys(y))

# Average results (if some have more runs).
for e in examples
    for k in keys(y[e])
        (k == "cnt") && continue
        y[e][k] /= y[e]["cnt"]
    end
end

# Average results (over all runs [already included in baseline], then down to "per iteration"), normalize to baseline.
e_0 = (60, 0, true)
baseline_iter = y[e_0]["iter"] / 100.0
baseline = (y[e_0]["main"] + y[e_0]["sub_p"]) / 100.0
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
        :yaxis_type => "log",
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
            y = 1.0,
            bordercolor = "black",
            borderwidth = 1,
            xanchor = "right",
            yanchor = "top",
            orientation = "v",
        ),
        :plot_bgcolor => "white",
        :paper_bgcolor => "white",
        :font => PlotlyJS.attr(; family = "Arial, sans-serif", size = 12, color = "black"),
        :margin => PlotlyJS.attr(; l = 60, r = 50, b = 65, t = 40),
    )

    return plot(traces, Layout(; kwlay..., layout...))
end

names = Dict(e => e[2] == 0 ? "never" : "$(e[2])" for e in examples)

for alg in ["ipm", "simplex"]
    exs = sort([e for e in examples if (e[3] ⊻ (alg == "simplex"))]; by = e -> e[2], rev = false)
    yn = [names[e] for e in exs]

    traces = Vector{PlotlyJS.GenericTrace}()

    zv = zeros(length(exs) - 1)
    for is_baseline in [true, false]
        alg == "simplex" && is_baseline && continue
        ib = is_baseline
        push!(
            traces,
            bar(;
                y = ib ? vcat(y[e_0]["iter"], zv) :
                    alg == "ipm" ? vcat(0, [y[e]["iter"] for e in exs if e != e_0]) :
                    [y[e]["iter"] for e in exs if e != e_0],
                x = yn,
                marker = !ib ? PlotlyJS.attr(; color = "#0f48aa") :
                         PlotlyJS.attr(;
                    color = "#0f48aa",
                    pattern_fillmode = "overlay",
                    pattern_fgcolor = "#ffffff",
                    pattern_bgcolor = "#0f48aa",
                    pattern_shape = "x",
                ),
                orientation = "v",
                name = "iterations",
                offsetgroup = 1,
                showlegend = !ib,
                # legendgrouptitle = PlotlyJS.attr(; text = "time"),
                # legendgroup = "iterations",
            ),
        )
        push!(
            traces,
            bar(;
                y = ib ? vcat(y[e_0]["main"] + y[e_0]["sub_s"], zv) :
                    alg == "ipm" ? vcat(0, [y[e]["main"] + y[e]["sub_s"] for e in exs if e != e_0]) :
                    [y[e]["main"] + y[e]["sub_s"] for e in exs if e != e_0],
                x = yn,
                marker = !ib ? PlotlyJS.attr(; color = "#458a0077") :
                         PlotlyJS.attr(;
                    color = "#458a0077",
                    pattern_fillmode = "overlay",
                    pattern_fgcolor = "#ffffff",
                    pattern_bgcolor = "#458a0077",
                    pattern_shape = "x",
                ),
                orientation = "v",
                name = "sub (serial)",
                offsetgroup = 2,
                legendgrouptitle = PlotlyJS.attr(; text = "model time"),
                legendgroup = "time",
                showlegend = !ib,
            ),
        )
        push!(
            traces,
            bar(;
                y = ib ? vcat(y[e_0]["main"] + y[e_0]["sub_p"], zv) :
                    alg == "ipm" ? vcat(0, [y[e]["main"] + y[e]["sub_p"] for e in exs if e != e_0]) :
                    [y[e]["main"] + y[e]["sub_p"] for e in exs if e != e_0],
                x = yn,
                orientation = "v",
                name = "sub (parallel)",
                offsetgroup = 2,
                marker = !ib ? PlotlyJS.attr(; color = "#458a00") :
                         PlotlyJS.attr(;
                    color = "#458a00",
                    pattern_fillmode = "overlay",
                    pattern_fgcolor = "#ffffff",
                    pattern_bgcolor = "#458a00",
                    pattern_shape = "x",
                ),
                legendgrouptitle = PlotlyJS.attr(; text = "model time"),
                legendgroup = "time",
                showlegend = !ib,
            ),
        )
        push!(
            traces,
            bar(;
                y = ib ? vcat(y[e_0]["main"], zv) :
                    alg == "ipm" ? vcat(0, [y[e]["main"] for e in exs if e != e_0]) :
                    [y[e]["main"] for e in exs if e != e_0],
                x = yn,
                marker = !ib ? PlotlyJS.attr(; color = "#b85c5c") :
                         PlotlyJS.attr(;
                    color = "#b85c5c",
                    pattern_fillmode = "overlay",
                    pattern_fgcolor = "#ffffff",
                    pattern_bgcolor = "#b85c5c",
                    pattern_shape = "x",
                ),
                orientation = "v",
                name = "main",
                offsetgroup = 2,
                legendgrouptitle = PlotlyJS.attr(; text = "time"),
                legendgroup = "time",
                showlegend = !ib,
            ),
        )
    end
    savefig(
        make_plot(
            traces,
            (
                barmode = "group",
                yaxis_title = "iterations / time (%)",
                xaxis_title = "drop non-binding cuts after n iterations",
                yaxis_range = [1, 3.2],
            ),
        ),
        joinpath(VIZ_DIR, "fig_$(alg).png");
        width = 600,
        height = 500,
    )
end

import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))

COLORS = ["#9c0000ff", "#009c00ff", "#00009cff"]
COLORS_T = ["#9c000077", "#009c0077", "#00009c77"]

PARALLELIZATION = 16
y = Dict()
y_keys = ["iter", "tmain_opt", "tmain_aux", "tsub_opt_ser", "tsub_aux_ser", "tsub_opt_par", "tsub_aux_par"]

files = []
for r in RUNS
    fh, ft = readdir(joinpath(RUN_DIR, r); join = true)
    i = parse(Int64, rsplit(rsplit(fh, "."; limit = 2)[1], "_"; limit = 2)[2])
    push!(files, (i, fh, ft))
end

for f in files
    n = f[1]
    haskey(y, n) || (y[n] = Dict{String, Any}(k => [] for k in y_keys))
    iter = JSON3.read(f[2]; allow_inf = true)[end]["k"] + 1
    timings = JSON3.read(f[3]; allow_inf = true)

    tim = timings[:inner_timers][:main][:inner_timers]
    tis = timings[:inner_timers][:sub][:inner_timers]

    tmain = tim[:optimize][:time_ns]
    tmain_aux = sum(v[:time_ns] for (k, v) in tim if k != :optimize)
    tsub = Dict(k => v[:inner_timers][:optimize][:time_ns] for (k, v) in tis)
    tsub_aux = Dict(sk => sum(v[:time_ns] for (k, v) in sv[:inner_timers] if k != :optimize) for (sk, sv) in tis)

    push!(y[n]["iter"], Float64(iter))
    push!(y[n]["tmain_opt"], tmain)
    push!(y[n]["tmain_aux"], tmain_aux)

    push!(y[n]["tsub_opt_ser"], sum(values(tsub)))
    push!(y[n]["tsub_aux_ser"], sum(values(tsub_aux)))

    par = [Dict("opt" => 0, "aux" => 0) for _ in 1:PARALLELIZATION]
    pi = 1
    for (_, i) in sort([(tsub[k] + tsub_aux[k], k) for k in keys(tsub)]; rev = true)
        par[pi]["opt"] += tsub[i]
        par[pi]["aux"] += tsub_aux[i]
        pi = pi % PARALLELIZATION + 1
    end
    push!(y[n]["tsub_opt_par"], par[1]["opt"])
    push!(y[n]["tsub_aux_par"], par[1]["aux"])
end

@info "EXPERIMENT $(EXPERIMENT_NR)" nof_experiments = length(y) avg_runs = length(files) / length(y)

# Average results and convert to seconds.
for n in keys(y)
    for k in keys(y[n])
        y[n][k] = mean(y[n][k]) / 1e9
    end
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

# x = string.(ts)
ts = sort(collect(keys(y)))
@assert all(ts .== [1, 4, 12, 24, 60, 365])
x = ["1Y", "1Q", "1M", "15D", "6D", "1D"]

# Plot iterations.
base = y[1]["iter"] / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["iter"] / base for n in ts], name = "main optimizer", marker_color = COLORS[3]))
savefig(
    make_plot(traces, (xaxis_title = "time-period of a single sub-model", yaxis_title = "iterations (%)")),
    joinpath(VIZ_DIR, "iterations.png");
    width = 400,
    height = 500,
)

# Plot time (serial).
base = sum(y[1][i] for i in ["tmain_opt", "tmain_aux", "tsub_opt_ser", "tsub_aux_ser"]) / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["tmain_opt"] / base for n in ts], name = "main (solve)", marker_color = COLORS[1]))
push!(
    traces,
    bar(; x = x, y = [y[n]["tmain_aux"] / base for n in ts], name = "main (overhead)", marker_color = COLORS_T[1]),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_opt_ser"] / base for n in ts], name = "sub (solve)", marker_color = COLORS[2]),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_aux_ser"] / base for n in ts], name = "sub (overhead)", marker_color = COLORS_T[2]),
)
savefig(
    make_plot(traces, (barmode = "stack", xaxis_title = "time-period of a single sub-model", yaxis_title = "time (%)")),
    joinpath(VIZ_DIR, "time_serial.png");
    width = 400,
    height = 500,
)

# Plot time (parallel).
base = sum(y[1][i] for i in ["tmain_opt", "tmain_aux", "tsub_opt_par", "tsub_aux_par"]) / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["tmain_opt"] / base for n in ts], name = "main (solve)", marker_color = COLORS[1]))
push!(
    traces,
    bar(; x = x, y = [y[n]["tmain_aux"] / base for n in ts], name = "main (overhead)", marker_color = COLORS_T[1]),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_opt_par"] / base for n in ts], name = "sub (solve)", marker_color = COLORS[2]),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_aux_par"] / base for n in ts], name = "sub (overhead)", marker_color = COLORS_T[2]),
)
savefig(
    make_plot(traces, (barmode = "stack", xaxis_title = "time-period of a single sub-model", yaxis_title = "time (%)")),
    joinpath(VIZ_DIR, "time_parallel.png");
    width = 400,
    height = 500,
)

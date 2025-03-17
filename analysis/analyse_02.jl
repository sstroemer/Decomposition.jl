import JSON3
using Statistics: mean
using PlotlyJS: PlotlyJS, plot, Layout, savefig, bar

EXPERIMENT_NR = split(split(basename(@__FILE__), ".")[1], "_")[2]
RESULT_DIR = normpath(@__DIR__, "..", "experiments", "out")
RUN_DIR = normpath(RESULT_DIR, only(filter(it -> startswith(it, "$(EXPERIMENT_NR)"), readdir(RESULT_DIR))))
RUNS = filter(x -> isdir(joinpath(RUN_DIR, x)), readdir(RUN_DIR))
VIZ_DIR = mkpath(replace(RUN_DIR, "experiments" => "analysis"))

PARALLELIZATION = 16
ts = [1, 4, 12, 24, 60, 365]
y = Dict(
    n => Dict{String, Any}(
        "iter" => [],
        "tmain_opt" => [],
        "tmain_aux" => [],
        "tsub_opt_ser" => [],
        "tsub_aux_ser" => [],
        "tsub_opt_par" => [],
        "tsub_aux_par" => [],
    ) for n in ts
)

# Extract results.
for r in RUNS
    dir = joinpath(RUN_DIR, r)

    for n in ts
        iter = JSON3.read(joinpath(dir, "history_$(n).json"))[end]["k"] + 1
        timings = JSON3.read(joinpath(dir, "timer_$(n).json"))

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
end

# Average results and convert to seconds.
for n in ts
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

x = string.(ts)

# Plot iterations.
base = y[1]["iter"] / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["iter"] / base for n in ts], name = "main optimizer", marker_color = "blue"))
savefig(
    make_plot(traces, (xaxis_title = "number of temporal splits", yaxis_title = "iterations (%)")),
    joinpath(VIZ_DIR, "iterations.png");
    width = 400,
    height = 500,
)

# Plot time (serial).
base = sum(y[1][i] for i in ["tmain_opt", "tmain_aux", "tsub_opt_ser", "tsub_aux_ser"]) / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["tmain_opt"] / base for n in ts], name = "main (solve)", marker_color = "#ff4800"))
push!(
    traces,
    bar(; x = x, y = [y[n]["tmain_aux"] / base for n in ts], name = "main (overhead)", marker_color = "#ffa480"),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_opt_ser"] / base for n in ts], name = "sub (solve)", marker_color = "#0026ff"),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_aux_ser"] / base for n in ts], name = "sub (overhead)", marker_color = "#8093ff"),
)
savefig(
    make_plot(traces, (barmode = "stack", xaxis_title = "number of temporal splits", yaxis_title = "time (%)")),
    joinpath(VIZ_DIR, "time_serial.png");
    width = 400,
    height = 500,
)

# Plot time (parallel).
base = sum(y[1][i] for i in ["tmain_opt", "tmain_aux", "tsub_opt_par", "tsub_aux_par"]) / 100.0
traces = Vector{PlotlyJS.GenericTrace}()
push!(traces, bar(; x = x, y = [y[n]["tmain_opt"] / base for n in ts], name = "main (solve)", marker_color = "#ff4800"))
push!(
    traces,
    bar(; x = x, y = [y[n]["tmain_aux"] / base for n in ts], name = "main (overhead)", marker_color = "#ffa480"),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_opt_par"] / base for n in ts], name = "sub (solve)", marker_color = "#0026ff"),
)
push!(
    traces,
    bar(; x = x, y = [y[n]["tsub_aux_par"] / base for n in ts], name = "sub (overhead)", marker_color = "#8093ff"),
)
savefig(
    make_plot(traces, (barmode = "stack", xaxis_title = "number of temporal splits", yaxis_title = "time (%)")),
    joinpath(VIZ_DIR, "time_parallel.png");
    width = 400,
    height = 500,
)

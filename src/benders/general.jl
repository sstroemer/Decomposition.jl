
abs_gap(x::Float64, y::Float64) = abs(x - y)
function rel_gap(x::Float64, y::Float64)
    tol = sqrt(eps(Float64))
    lower = min(x, y)
    upper = max(x, y)
    gap = upper - lower

    # This is similar to: https://www.gurobi.com/documentation/current/refman/mipgap2.html
    # Note: Behaviour for `upper == 0` may not be as expected.

    (isnan(gap) || !isfinite(gap)) && return +Inf
    isapprox(upper, 0.0; atol=tol) && return isapprox(lower, 0.0; atol=tol) ? 0.0 : +Inf
    return abs(gap / upper)
end

main(model::Benders.DecomposedModel) = model.models[1]
sub(model::Benders.DecomposedModel; index::Int) = model.models[1 + index]
subs(model::Benders.DecomposedModel) = model.models[2:end]

current_iteration(model::Benders.DecomposedModel) = length(model.info[:history])
best_upper_bound(model::Benders.DecomposedModel) = minimum(it[:upper_bound] for it in model.info[:history]; init=+Inf)
best_lower_bound(model::Benders.DecomposedModel) = maximum(it[:lower_bound] for it in model.info[:history]; init=-Inf)
best_gap_abs(model::Benders.DecomposedModel) = abs_gap(best_lower_bound(model), best_upper_bound(model))
best_gap_rel(model::Benders.DecomposedModel) = rel_gap(best_lower_bound(model), best_upper_bound(model))

total_wall_time(model::Benders.DecomposedModel) = sum(it[:time][:wall] for it in model.info[:history]; init=0)
total_cpu_time(model::Benders.DecomposedModel) = sum(it[:time][:cpu] for it in model.info[:history]; init=0)

get_main_time(model::Benders.DecomposedModel) = TimerOutputs.time(model.timer["main"])
get_sub_time(model::Benders.DecomposedModel; index::Int) = TimerOutputs.time(model.timer["sub"]["[$(index)]"])

function next_iteration!(model::Benders.DecomposedModel, added_cuts; verbose::Bool=true, assumed_pcores::Int = 16)
    nof_feas_cuts, nof_opt_cuts = added_cuts
    iter = current_iteration(model)
    curr_time = time_ns()

    # Calculate the lower bound.
    lb = model.info[:results][:main][:obj_f_val] # TODO: should be `res_main_obj`? see: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592

    # Calculate the upper bound.
    subs_obj = [it[:obj] for it in model.info[:results][:subs]]
    ub = any(ismissing, subs_obj) ? Inf : sum(subs_obj)
    ub += model.info[:results][:main][:obj_f_val] - sum(model.info[:results][:main][:θ])
    # TODO: Account here for stuff like "ObjectiveInCuts"

    # Find all cuts that were added in this iteration.
    new_cuts = OrderedDict(
        :feasibility => model.cuts[:feasibility][(end-nof_feas_cuts+1):end],
        :optimality => model.cuts[:optimality][(end-nof_opt_cuts+1):end],
    )

    # Estimate "batched" parallel processing of sub-model timings.
    # TODO: Reimplement this for TimerOutputs
    time_est = get_main_time(model) + sum(get_sub_time(model, index=i) for i in 1:length(subs(model))) - total_wall_time(model)
    # _t = est_Δt_wall[:subs]
    # _batches = [
    #     _t[((i-1) * assumed_pcores + 1):min(i * assumed_pcores, end)]
    #     for i in 1:ceil(Int, length(_t) / assumed_pcores)
    # ]
    # subs_wall_time = sum(maximum.(_batches); init=Inf)

    # Prepare and add the history entry.
    entry = OrderedDict(
        # TODO: include "quality criteria":
        #       - whether all sub-models are: (1) feasible, (2) w/o the use of artificial slacks
        #       - conditioning, ... of main and sub-models
        #       - ...
        :iteration => iter,
        :time => OrderedDict(
            :timestamp => curr_time,
            :wall => time_est, # TODO
            :cpu => time_est,  # TODO
            :estimate => time_est,  # TODO
        ),
        :lower_bound => lb,
        :upper_bound => ub,
        :gap_abs => abs_gap(lb, ub),
        :gap_rel => rel_gap(lb, ub),
        :added_cuts => OrderedDict(:feasibility => nof_feas_cuts, :optimality => nof_opt_cuts),
        :added_cuts_con => new_cuts,
    )
    push!(model.info[:history], entry)

    # Printing.
    if verbose
        if iter == 0
            # Print motd-like header.
            motd = [
                "╭──────────────────────────────────────────────────────────────────────────────────────────────────────────────╮",
                "│ >> Decomposition.jl <<                                              [version::0.1.0  //  algorithm::benders] │",
                "├────────┬───────────────────────────┬───────────────────────────┬───────────────────────────┬─────────────────┤",
                "│        │      objective bound      │      current best gap     │    est. execution time    │   added cuts    │",
                "├────────┼─────────────┬─────────────┼───────────────────────────┤─────────────┬─────────────┤────────┬────────┤",
                "│   iter │       lower │       upper │    absolute │    relative │    wall [s] │     cpu [s] │  feas. │   opt. │",
                "├────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼────────┼────────┤",  
            ]
            append!(model.log, motd)
            println.(motd)
        end

        push!(
            model.log,
            _print_iteration(
                iter,
                best_lower_bound(model),
                best_upper_bound(model),
                best_gap_abs(model),
                best_gap_rel(model),
                Printf.@sprintf("%11.2f", total_wall_time(model) / 1e9),
                Printf.@sprintf("%11.2f", total_cpu_time(model) / 1e9),
                length(model.cuts[:feasibility]),
                length(model.cuts[:optimality]),
            )
        )
    end

    return nothing
end

function save(model::Benders.DecomposedModel)
    sio = IOBuffer()
    JSON3.pretty(
        sio,
        OrderedDict(
            "meta" => OrderedDict(
                "version" => "0.1.0",
                "algorithm" => "benders",
                "name" => model.name,
                "created" => model.info[:stats][:created],
                "started" => model.info[:stats][:started],
            ),
            "inputs" => OrderedDict(
                "timesteps" => model.T,
                "splits" => model.nof_temporal_splits,
            ),
            "attributes" => showtostr.(model.attributes),
            "models" => OrderedDict(
                "main" => _showtostr(main(model)),
                "subs" => _showtostr.(subs(model)),
            ),
            "cuts" => OrderedDict(
                "feasibility" => length(model.cuts[:feasibility]),
                "optimality" => length(model.cuts[:optimality]),
            ),
            "history" => filter.(k -> (k.first != :added_cuts_con), model.info[:history]),
            "log" => join(model.log, "\n"),
        ),
        JSON3.AlignmentContext(alignment=:Left, indent=2);
        allow_inf=true
    )
    
    json_str = String(take!(sio))
    short_hash = bytes2hex(SHA.sha1(json_str))[1:7]
    filename = "$(model.name)_$(short_hash).djl.json"

    open(normpath(mkpath("out"), filename), "w") do f
        write(f, json_str)
        @info "Saved model" filename
    end

    return nothing
end

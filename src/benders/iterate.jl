function iterate!(model::Benders.DecomposedModel; nthreads::Int = -1)   
    @timeit model.timer "main" begin
        @timeit model.timer "optimize" execute!(model, Benders.Query.SolveMain())
        @timeit model.timer "extract solution" execute!(model, Benders.Query.ExtractResultsMain())
        nothing
    end

    # TODO: we could abort here if the gap is small enough (saving one final sub-model iteration) ?

    # TODO: reactivate and refactor sub-batching

    # sub_batching = nothing
    # if nthreads > 1
    #     @timeit model.timer "main" begin
    #             sub_batching = @timeit model.timer "orchestrate" Benders.allocate_sub_models(model, nthreads)
    #     end
    # end

    @timeit model.timer "sub" begin
        # TODO: also remove this and reactivate below
        for i in 1:(length(model.models) - 1)
            @timeit model.timer "[$i]" begin
                execute!(model, Benders.Query.SolveSub(index=i))
                execute!(model, Benders.Query.ExtractResultsSub(index=i))
            end
        end

        # if nthreads < 1
        #     for i in 1:(length(model.models) - 1)
        #         @timeit model.timer "[$i]" begin
        #             Benders.solve_sub(model; index=i)
        #             Benders.extract_sub(model; index=i)
        #         end
        #     end
        # else
        #     Threads.@threads for batch in sub_batching
        #         for i in batch
        #             @timeit model.timer "[$i]" begin
        #                 Benders.solve_sub(model; index=i)
        #                 Benders.extract_sub(model; index=i)
        #             end
        #         end
        #     end
        # end
    end

    added_cuts = @timeit model.timer "main" begin
        new_cuts = @timeit model.timer "cuts (generate)" generate_cuts(model)
        @timeit model.timer "cuts (preprocess)" preprocess_cuts!(model, new_cuts)
        added_cuts = @timeit model.timer "cuts (add)" add_cuts(model, new_cuts)
        @timeit model.timer "cuts (postprocess)" postprocess_cuts!(model)

        added_cuts
    end
    
    # Pass added cuts to `next_iteration!`.
    next_iteration!(model, added_cuts)

    return check_termination(model)
end

function next_iteration!(model::Benders.DecomposedModel, added_cuts; verbose::Bool=true, assumed_pcores::Int = 16)
    nof_feas_cuts, nof_opt_cuts = added_cuts
    iter = current_iteration(model)
    curr_time = time_ns()

    # TODO: move the 2 steps afterwards into a separate function "calc_global_bounds"

    # Calculate the lower bound.
    lb = coalesce(model.info[:results][:main][:obj_lb], -Inf)

    # Calculate the upper bound.
    # TODO: after refactoring "extract" for subs, the below should be the same for all cut types
    if has_attribute_type(model, Benders.CutTypeMISFSZ)
        ub = 0.0
        ub += model.info[:results][:main][:obj_base]
        ub += coalesce(sum(res[:obj_val_primal] for res in model.info[:results][:subs]), +Inf)
    else
        subs_obj = [it[:obj] for it in model.info[:results][:subs]]
        ub = any(ismissing, subs_obj) ? Inf : sum(subs_obj)
        ub += model.info[:results][:main][:obj_base]
    end

    # Find all cuts that were added in this iteration.
    if has_attribute_type(model, Benders.CutTypeMISFSZ)
        new_cuts = OrderedDict(
            :feasibility => Benders.AbstractGeneralCut[],
            :optimality => Benders.AbstractGeneralCut[],
            :misfsz => model.cuts[:misfsz][(end-nof_feas_cuts-nof_opt_cuts+1):end],
        )
    else
        new_cuts = OrderedDict(
            :feasibility => model.cuts[:feasibility][(end-nof_feas_cuts+1):end],
            :optimality => model.cuts[:optimality][(end-nof_opt_cuts+1):end],
            :misfsz => Benders.AbstractGeneralCut[],
        )
    end

    # Estimate "batched" parallel processing of sub-model timings.
    # TODO: Reimplement this for TimerOutputs
    time_est = (
        get_main_time(model) +
        sum(get_sub_time(model, index=i) for i in 1:length(Benders.subs(model))) -
        total_wall_time(model)
    )
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
        :gap_abs => gap_abs(lb, ub),
        :gap_rel => gap_rel(lb, ub),
        :added_cuts => has_attribute_type(model, Benders.CutTypeMISFSZ) ? OrderedDict(:feasibility => 0, :optimality => 0, :misfsz => nof_feas_cuts+nof_opt_cuts) : OrderedDict(:feasibility => nof_feas_cuts, :optimality => nof_opt_cuts, :misfsz => 0),
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
                length(model.cuts[:feasibility]) + count(c -> c isa Benders.MISFSZFeasibilityCut, model.cuts[:misfsz]),
                length(model.cuts[:optimality]) + count(c -> c isa Benders.MISFSZOptimalityCut, model.cuts[:misfsz]),
            )
        )
    end

    # Clean "results" for next iteration.
    model.info[:results] = OrderedDict{Symbol, Any}(
        :main => OrderedDict{Symbol, Any}(),
        :subs => [OrderedDict{Symbol, Any}() for _ in Benders.subs(model)],
    )

    return nothing
end

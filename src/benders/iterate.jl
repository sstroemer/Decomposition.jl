function iterate!(model::Benders.DecomposedModel; nthreads::Int = -1)   
    @timeit model.timer "main" begin
        @timeit model.timer "optimize" solve_main(model)
        @timeit model.timer "extract solution" extract_main(model)
    end

    # TODO: we could abort here if the gap is small enough (saving one final sub-model iteration) ?

    sub_batching = nothing
    if nthreads > 1
        @timeit model.timer "main" begin
                sub_batching = @timeit model.timer "orchestrate" allocate_sub_models(model, nthreads)
        end
    end

    @timeit model.timer "sub" begin
        if nthreads < 1
            for i in 1:(length(model.models) - 1)
                @timeit model.timer "[$i]" begin
                    solve_sub(model; index=i)
                    extract_sub(model; index=i)
                end
            end
        else
            Threads.@threads for batch in sub_batching
                for i in batch
                    @timeit model.timer "[$i]" begin
                        solve_sub(model; index=i)
                        extract_sub(model; index=i)
                    end
                end
            end
        end
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

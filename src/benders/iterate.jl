function iterate!(model::Benders.DecomposedModel)   
    @timeit model.timer "main" begin
        @timeit model.timer "optimize" solve_main(model)
        @timeit model.timer "extract solution" extract_main(model)
    end

    # TODO: we could abort here if the gap is small enough (saving one final sub-model iteration) ?

    @timeit model.timer "sub" begin
        for i in 1:(length(model.models) - 1)
            @timeit model.timer "[$i]" begin
                solve_sub(model; index=i)
                extract_sub(model; index=i)
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

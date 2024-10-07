

function estimate_sub_model_time(model::Benders.DecomposedModel)
    function _estimate(model::Benders.DecomposedModel; index::Int)
        if haskey(model.timer, "sub") && haskey(model.timer["sub"], "[$(index)]")
            return model.timer["sub"]["[$(index)]"].accumulated_data.time / model.timer["sub"]["[$(index)]"].accumulated_data.ncalls
        else
            return 1.0  # TODO: base that on number of vars/cons/nnzs, etc.
        end
    end

    return sort(collect(Dict(i => _estimate(model; index=i) / 1e9 for i in 1:length(Benders.subs(model)))); by=x->x[2], rev=true)
end

function allocate_sub_models(model::Benders.DecomposedModel, n::Int)
    if n == 1
        return collect(1:length(Benders.subs(model)))
    end

    batches = [Int64[] for _ in 1:n]
    batch_times = [0.0 for _ in 1:n]

    estimated_time = estimate_sub_model_time(model)

    for (i, t) in estimated_time
        batch_idx = argmin(batch_times)
        push!(batches[batch_idx], i)
        batch_times[batch_idx] += t
    end

    return batches
end

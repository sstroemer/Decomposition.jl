abstract type DecompositionAttribute end
# TODO: Keep a list of all "attributes" (which are just sub-types that describe modifications, and may contain more data (e.g., the sub-model index they were applied to))

@kwdef struct DecomposedModel1
    monolithic::JuMP.Model

    lpmd::JuMP.LPMatrixData

    idx_model_vars::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]
    idx_model_cons::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]


    # TODO: allow storing "all_variables" somewhere for each model (and for main: before adding θ!)
    models::Vector{JuMP.Model} = Vector{JuMP.Model}[]
    reference_maps = Vector{Any}()

    annotations = Dict{Any, Any}(:variables => Dict{Any, Any}(), :constraints => Dict{Any, Any}())

    decomposition_maps = Dict{Any, Any}()

    stats = Dict{Symbol, Any}(
        :iteration => 0,
        :lower_bound => -Inf,
        :upper_bound => Inf,
        :gap_rel => Inf,
        :gap_abs => Inf,
    )     # TODO: transform into a struct, that tracks stats for each iteration (using a "inc_iter" function)

    cuts = Vector{Any}()  # TODO: track which cut is created by which iteration (inside the stats struct)
end
DecomposedModel = DecomposedModel1

attach_solver(model::JuMP.Model, solver) = set_optimizer(model, solver)  # TODO: use this to attach "BD" (or others)
solve!(model::JuMP.Model) = optimize!(model)

function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64})
    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too
    # TODO: transform single variable constraints into bounds

    model = direct_model(HiGHS.Optimizer())

    # Create variables, and set bounds.
    @variable(model, x[i = idx_v])
    for i in eachindex(x)
        isfinite(lpmd.x_lower[i]) && set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && set_upper_bound(x[i], lpmd.x_upper[i])
    end

    # Create constraints.
    for i in idx_c
        if lpmd.b_lower[i] == lpmd.b_upper[i]
            # If `lb == ub`, then we know both have to be finite.
            @constraint(model, lpmd.A[i, idx_v]' * x.data == lpmd.b_lower[i])
            continue
        end

        isfinite(lpmd.b_lower[i]) && @constraint(model, lpmd.A[i, idx_v]' * x.data >= lpmd.b_lower[i])
        isfinite(lpmd.b_upper[i]) && @constraint(model, lpmd.A[i, idx_v]' * x.data <= lpmd.b_upper[i])
    end

    return model
end

abstract type DecompositionAttribute end
abstract type DecompositionQuery end

@kwdef struct DecomposedModel6
    monolithic::JuMP.Model
    lpmd::JuMP.LPMatrixData

    T::Int64
    nof_temporal_splits::Int64

    # TODO: Store these as (sorted) vector (is that even needed?) and as set (for faster lookup)
    idx_model_vars::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]
    idx_model_cons::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]


    # TODO: allow storing "all_variables" somewhere for each model (and for main: before adding θ!)
    models::Vector{JuMP.Model} = Vector{JuMP.Model}[]
    reference_maps = Vector{Any}()

    annotations = Dict{Any, Any}(:variables => Dict{Any, Any}(), :constraints => Dict{Any, Any}())

    decomposition_maps = Dict{Any, Any}()

    info = Dict{Symbol, Any}(
        :stats => Dict(
            :created => time_ns(),
            :started => missing,
        ),
        :history => Vector{Dict{Symbol, Any}}(),
        :results => Dict{Symbol, Any}(
            :main => Dict{Symbol, Any}(),
            :subs => Vector{Dict}()
        )
    )

    attributes::Vector{DecompositionAttribute} = DecompositionAttribute[]

    cuts = Dict{Symbol, Vector{ConstraintRef}}(
        :feasibility => ConstraintRef[],
        :optimality => ConstraintRef[],
    )  # TODO: track which cut is created by which iteration (inside the stats struct)
end
DecomposedModel = DecomposedModel6

_abs_gap(x::Float64, y::Float64) = abs(x - y)
function _rel_gap(x::Float64, y::Float64)
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

current_iteration(model::DecomposedModel) = length(model.info[:history])
best_upper_bound(model::DecomposedModel) = minimum(it[:upper_bound] for it in model.info[:history]; init=+Inf)
best_lower_bound(model::DecomposedModel) = maximum(it[:lower_bound] for it in model.info[:history]; init=-Inf)
best_gap_abs(model::DecomposedModel) = _abs_gap(best_lower_bound(model), best_upper_bound(model))
best_gap_rel(model::DecomposedModel) = _rel_gap(best_lower_bound(model), best_upper_bound(model))

total_wall_time(model::DecomposedModel) = sum(it[:time][:wall] for it in model.info[:history])
total_cpu_time(model::DecomposedModel) = sum(it[:time][:cpu] for it in model.info[:history])

function next_iteration!(model::DecomposedModel, added_cuts, est_Δt_wall; verbose::Bool=true, assumed_pcores::Int = 16)
    nof_feas_cuts, nof_opt_cuts = added_cuts
    iter = current_iteration(model)
    curr_time = time_ns()

    # Calculate the lower bound.
    lb = model.info[:results][:main][:obj_f_val] # TODO: should be `res_main_obj`? see: https://discourse.julialang.org/t/detecting-problems-with-numerically-challenging-models/118592

    # Calculate the upper bound.
    subs_obj = [it[:obj] for it in model.info[:results][:subs]]
    ub = any(ismissing, subs_obj) ? Inf : sum(subs_obj)
    if bd_has_attribute_type(model, BD_MainObjectiveCon)
        ub += model.info[:results][:main][:obj_exp]
    else
        ub += model.info[:results][:main][:obj_f_val] - sum(model.info[:results][:main][:θ])
    end

    # Find all cuts that were added in this iteration.
    new_cuts = Dict(
        :feasibility => model.cuts[:feasibility][(end-nof_feas_cuts+1):end],
        :optimality => model.cuts[:optimality][(end-nof_opt_cuts+1):end],
    )

    # Estimate "batched" parallel processing of sub-model timings.
    _t = est_Δt_wall[:subs]
    _batches = [
        _t[((i-1) * assumed_pcores + 1):min(i * assumed_pcores, end)]
        for i in 1:ceil(Int, length(_t) / assumed_pcores)
    ]
    subs_wall_time = sum(maximum.(_batches))

    # Prepare and add the history entry.
    entry = Dict(
        :iteration => iter,
        :time => Dict(
            :timestamp => curr_time,
            :wall => est_Δt_wall[:main] + subs_wall_time + est_Δt_wall[:aux],
            :cpu => curr_time - (iter == 0 ? model.info[:stats][:started] : model.info[:history][end][:time][:timestamp]),
            :estimate => est_Δt_wall,
        ),
        :lower_bound => lb,
        :upper_bound => ub,
        :gap_abs => _abs_gap(lb, ub),
        :gap_rel => _rel_gap(lb, ub),
        :added_cuts => new_cuts,
    )
    push!(model.info[:history], entry)

    # Printing.
    if verbose
        if iter == 0
            # Print motd-like header.
            println("╭──────────────────────────────────────────────────────────────────────────────────────────────────────────────╮")
            println("│ >> Decomposition.jl <<                                              [version::0.1.0  //  algorithm::benders] │")
            println("├────────┬───────────────────────────┬───────────────────────────┬───────────────────────────┬─────────────────┤")
            println("│        │      objective bound      │      current best gap     │    est. execution time    │   added cuts    │")
            println("├────────┼─────────────┬─────────────┼───────────────────────────┤─────────────┬─────────────┤────────┬────────┤")
            println("│   iter │       lower │       upper │    absolute │    relative │    wall [s] │     cpu [s] │  feas. │   opt. │")
            println("├────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼─────────────┼────────┼────────┤")
        end

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
    end

    return nothing
end

function check_termination(model::DecomposedModel)
    # TODO: implement
    return false
end

function modify(::DecomposedModel, ::DecompositionAttribute)
    @error "Not implemented"
end

attach_solver(model::JuMP.Model, solver) = set_optimizer(model, solver)  # TODO: use this to attach "BD" (or others)
solve!(model::JuMP.Model) = optimize!(model)

function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64})
    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too
    # TODO: transform single variable constraints into bounds

    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # model = direct_model(HiGHS.Optimizer())
    set_silent(model)

    # Create variables, and set bounds.
    @variable(model, x[i = idx_v])
    for i in idx_v
        isfinite(lpmd.x_lower[i]) && set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && set_upper_bound(x[i], lpmd.x_upper[i])
    end

    single_var_con = Set{Int64}()
    srv = sort(lpmd.A.rowval)
    for i in 3:length(srv)
        srv[i] == srv[i-1] && continue
        srv[i-1] == srv[i-2] &&  continue
        push!(single_var_con, srv[i-1])      
    end

    nzA = lpmd.A .!= 0
    fnzc = nzA * collect(1:size(nzA, 2))

    # Create constraints.
    add_ge = Vector{Int64}()
    add_lt = Vector{Int64}()
    add_eq = Vector{Int64}()
    for i in idx_c
        # Check for "single variable constraints", which can be transformed into bounds.
        if i in single_var_con
            _idx = fnzc[i]
            c = lpmd.A[i, _idx]
            
            if isfinite(lpmd.b_lower[i])
                new_lb = lpmd.b_lower[i] / c
                if !has_lower_bound(x[_idx]) || new_lb > lower_bound(x[_idx])
                    set_lower_bound(x[_idx], new_lb)
                end
            end

            if isfinite(lpmd.b_upper[i])
                new_ub = lpmd.b_upper[i] / c
                if !has_upper_bound(x[_idx]) || new_ub < upper_bound(x[_idx])
                    set_upper_bound(x[_idx], new_ub)
                end
            end

            continue
        end

        if lpmd.b_lower[i] == lpmd.b_upper[i]
            # If `lb == ub`, then we know both have to be finite.
            push!(add_eq, i)
            continue
        end

        isfinite(lpmd.b_lower[i]) && push!(add_ge, i)
        isfinite(lpmd.b_upper[i]) && push!(add_lt, i)
    end

    @constraint(model, lpmd.A[add_eq, idx_v] * x.data .== lpmd.b_lower[add_eq])
    @constraint(model, lpmd.A[add_ge, idx_v] * x.data .>= lpmd.b_lower[add_ge])
    @constraint(model, lpmd.A[add_lt, idx_v] * x.data .<= lpmd.b_upper[add_lt])

    return model
end

struct SOLVE_AlgorithmSimplex <: DecompositionAttribute
    model::Symbol   # :main, :sub
    type::Symbol    # :primal, :dual, ...
end

struct SOLVE_AlgorithmIPM <: DecompositionAttribute
    model::Symbol   # :main, :sub
    crossover::Bool

    SOLVE_AlgorithmIPM(model::Symbol) = new(model, false)
end

function modify(model::DecomposedModel, attribute::SOLVE_AlgorithmSimplex)
    models = attribute.model == :main ? [bd_main(model)] : bd_subs(model)
    for m in models
        solver = solver_name(m)
        if solver == "Gurobi"
            val = -1
            if attribute.type == :primal
                val = 0
            elseif attribute.type == :dual
                val = 1
            else
                @error "Simplex type `$(attribute.type)` not supported by solver `$(solver)`"
            end

            set_attribute(m, "Method", val)
        else
            @error "Setting `SOLVE_AlgorithmSimplex` is currently not supported for solver `$(solver)`"
        end
    end
end

function modify(model::DecomposedModel, attribute::SOLVE_AlgorithmIPM)
    models = attribute.model == :main ? [bd_main(model)] : bd_subs(model)
    for m in models
        solver = solver_name(m)
        if solver == "Gurobi"
            set_attribute(m, "Method", 2)
            set_attribute(m, "Crossover", attribute.crossover ? -1 : 0)
        elseif solver == "HiGHS"
            set_attribute(m, "solver", "ipm")
            set_attribute(m, "run_crossover", attribute.crossover ? "on" : "off")
        else
            @error "Setting `SOLVE_AlgorithmSimplex` is currently not supported for solver `$(solver)`"
        end
    end
end

"""
Helper function to get the graph structure.
"""
function create_adjacency_matrix(A::SparseArrays.SparseMatrixCSC, rm::Set{Int64})
    tA = SparseArrays.SparseMatrixCSC((A .!= 0)')
    rvs = SparseArrays.rowvals(tA)
    
    I = Int64[]
    J = Int64[]

    for i in 1:size(A, 1)
        nodes = rvs[SparseArrays.nzrange(tA, i)] # tA[:, i].nzind
        for j in eachindex(nodes)
            u = nodes[j]
            (u in rm) && continue
            for k in (j+1):length(nodes)
                v = nodes[k]
                (v in rm) && continue
                push!(I, u)
                push!(J, v)
            end
        end
    end

    return LinearAlgebra.Symmetric(SparseArrays.sparse(I, J, 1, size(A, 2), size(A, 2)))
end
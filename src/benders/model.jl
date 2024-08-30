@kwdef struct DecomposedModel1 <: AbstractDecomposedModel
    name::String

    monolithic::JuMP.Model
    lpmd::JuMP.LPMatrixData

    T::Int64
    nof_temporal_splits::Int64

    # TODO: Store these as (sorted) vector (is that even needed?) and as set (for faster lookup)
    vis::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]
    cis::Vector{Vector{Int64}} = Vector{Vector{Int64}}[]

    # TODO: Allow storing "all_variables" somewhere for each model (and for main: before adding θ!)
    models::Vector{JuMP.Model} = Vector{JuMP.Model}[]

    annotations = Dict{Any, Any}(:variables => Dict{Any, Any}(), :constraints => Dict{Any, Any}())

    attributes::Vector{AbstractDecompositionAttribute} = AbstractDecompositionAttribute[]
    _attribute_iteration_info::Vector{Int64} = Int64[]

    solutions = Dict{Symbol, Any}(
        :current => Dict{Symbol, Any}(
            :main => nothing
            :subs => nothing
        ),
    )

    # TODO: Move that into a struct
    info = OrderedDict{Symbol, Any}(
        :stats => OrderedDict(
            :created => time_ns(),
            :started => missing,
        ),
        :history => Vector{Dict{Symbol, Any}}(),
        :results => OrderedDict{Symbol, Any}(
            # TODO: merge this (that only holds "objective values") into `solutions`
            :main => OrderedDict{Symbol, Any}(),
            :subs => Vector{Dict}()
        )
    )

    cuts = OrderedDict{Symbol, Vector{JuMP.ConstraintRef}}(
        :feasibility => JuMP.ConstraintRef[],
        :optimality => JuMP.ConstraintRef[],
    )

    f_opt_main = () -> Gurobi.Optimizer(GRB_ENV)
    f_opt_sub = () -> Gurobi.Optimizer(GRB_ENV)

    log::Vector{String} = String[]
end
DecomposedModel = DecomposedModel1

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

current_iteration(model::DecomposedModel) = length(model.info[:history])
best_upper_bound(model::DecomposedModel) = minimum(it[:upper_bound] for it in model.info[:history]; init=+Inf)
best_lower_bound(model::DecomposedModel) = maximum(it[:lower_bound] for it in model.info[:history]; init=-Inf)
best_gap_abs(model::DecomposedModel) = abs_gap(best_lower_bound(model), best_upper_bound(model))
best_gap_rel(model::DecomposedModel) = rel_gap(best_lower_bound(model), best_upper_bound(model))

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
    ub += model.info[:results][:main][:obj_f_val] - sum(model.info[:results][:main][:θ])
    # TODO: Account here for stuff like "ObjectiveInCuts"

    # Find all cuts that were added in this iteration.
    new_cuts = OrderedDict(
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
    entry = OrderedDict(
        # TODO: include "quality criteria":
        #       - whether all sub-models are: (1) feasible, (2) w/o the use of artificial slacks
        #       - conditioning, ... of main and sub-models
        #       - ...
        :iteration => iter,
        :time => OrderedDict(
            :timestamp => curr_time,
            :wall => est_Δt_wall[:main] + subs_wall_time + est_Δt_wall[:aux],
            :cpu => curr_time - (iter == 0 ? model.info[:stats][:started] : model.info[:history][end][:time][:timestamp]),
            :estimate => est_Δt_wall,
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

function check_termination(model::DecomposedModel)
    # TODO: implement
    return false
end

function save(model::DecomposedModel)
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
                "main" => showtostr(main(model)),
                "subs" => showtostr.(subs(model)),
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

function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64}; optimizer)
    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too

    model = JuMP.direct_model(optimizer)
    # model = Model(() -> optimizer)
    JuMP.set_silent(model)

    # Create variables, and set bounds.
    JuMP.@variable(model, x[i = idx_v])
    for i in idx_v
        isfinite(lpmd.x_lower[i]) && JuMP.set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && JuMP.set_upper_bound(x[i], lpmd.x_upper[i])
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
                if !JuMP.has_lower_bound(x[_idx]) || new_lb > JuMP.lower_bound(x[_idx])
                    JuMP.set_lower_bound(x[_idx], new_lb)
                end
            end

            if isfinite(lpmd.b_upper[i])
                new_ub = lpmd.b_upper[i] / c
                if !JuMP.has_upper_bound(x[_idx]) || new_ub < JuMP.upper_bound(x[_idx])
                    JuMP.set_upper_bound(x[_idx], new_ub)
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

    JuMP.@constraint(model, lpmd.A[add_eq, idx_v] * x.data .== lpmd.b_lower[add_eq])
    JuMP.@constraint(model, lpmd.A[add_ge, idx_v] * x.data .>= lpmd.b_lower[add_ge])
    JuMP.@constraint(model, lpmd.A[add_lt, idx_v] * x.data .<= lpmd.b_upper[add_lt])

    return model
end

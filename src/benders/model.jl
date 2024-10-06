@kwdef struct DecomposedModel5 <: AbstractDecomposedModel
    jump_model::JuMP.Model
    T::Int64
    nof_temporal_splits::Int64

    name::String = pop!(jump_model.ext, :_model_name, "unnamed")

    # -------------

    f_opt::Base.Callable = () -> @error "No optimizer provided; either set `f_opt`, or `f_opt_main` and `f_opt_sub`"
    f_opt_main::Base.Callable = f_opt
    f_opt_sub::Base.Callable = f_opt

    # -------------

    lpmd::JuMP.LPMatrixData = JuMP.lp_matrix_data(jump_model)

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
            :main => nothing,
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
            :subs => Vector{OrderedDict{Symbol, Any}}(),
        )
    )

    cuts = OrderedDict{Symbol, Vector{AbstractCut}}(
        :feasibility => AbstractCut[],
        :optimality => AbstractCut[],
    )

    log::Vector{String} = String[]
    timer::TimerOutput = TimerOutput()
end
DecomposedModel = DecomposedModel5

Base.Broadcast.broadcastable(model::DecomposedModel) = Ref(model)
Base.show(io::IO, model::DecomposedModel) = print(io, "DecomposedModel [algorithm=Benders]: $(model.name)")

function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64}; optimizer, cache::Dict{Symbol, Any})
    # TODO: pull "debug", "verbosity", and other settings from the list of (MOI.RawOptimizationAttribute) attributes that were added to the model
    debug = false

    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too

    model = JuMP.direct_model(optimizer)
    # model = Model(() -> optimizer)
    JuMP.set_silent(model)
    JuMP.set_string_names_on_creation(model, debug)

    # Create variables, and set bounds.
    JuMP.@variable(model, x[i = idx_v])
    for i in idx_v
        isfinite(lpmd.x_lower[i]) && JuMP.set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && JuMP.set_upper_bound(x[i], lpmd.x_upper[i])
    end

    # Create constraints.
    add_ge = Vector{Int64}()
    add_lt = Vector{Int64}()
    add_eq = Vector{Int64}()
    for i in idx_c
        # Check for "single variable constraints", which can be transformed into bounds.
        if i in cache[:single_var_con]
            _idx = cache[:fnzc][i]
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

    # Instead of:
    # ```julia
    # JuMP.@constraint(model, (@view lpmd.A[add_eq, idx_v]) * x.data .== lpmd.b_lower[add_eq])
    # JuMP.@constraint(model, (@view lpmd.A[add_ge, idx_v]) * x.data .>= lpmd.b_lower[add_ge])
    # JuMP.@constraint(model, (@view lpmd.A[add_lt, idx_v]) * x.data .<= lpmd.b_upper[add_lt])
    # ```
    #
    # we optimize to:
    # ```julia
    # x_data_t = collect(x.data')
    # A_t = copy(lpmd.A[:, idx_v]')
    # 
    # isempty(add_eq) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_eq], 1.0)' .== lpmd.b_lower[add_eq])
    # isempty(add_ge) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_ge], 1.0)' .>= lpmd.b_lower[add_ge])
    # isempty(add_lt) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_lt], 1.0)' .<= lpmd.b_upper[add_lt])
    # ```
    # which again can be optimized by only transpose-copying once:

    x_data_t = collect(x.data')
    A_t = copy(lpmd.A[vcat(add_eq, add_ge, add_lt), idx_v]')  # TODO: this line takes > 50% of the total time
    
    idx_eq = 1:length(add_eq)
    idx_ge = length(add_eq)+1:length(add_eq)+length(add_ge)
    idx_lt = length(add_eq)+length(add_ge)+1:length(add_eq)+length(add_ge)+length(add_lt)   

    isempty(add_eq) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, idx_eq], 1.0)' .== lpmd.b_lower[add_eq])
    isempty(add_ge) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, idx_ge], 1.0)' .>= lpmd.b_lower[add_ge])
    isempty(add_lt) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, idx_lt], 1.0)' .<= lpmd.b_upper[add_lt])

    # Add default processes.
    model.ext[:processes] = Dict(
        :solve => [() -> JuMP.optimize!(model)],
        :extract => [],
    )

    return model
end

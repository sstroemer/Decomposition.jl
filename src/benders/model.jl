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
            :subs => Vector{Dict}()
        )
    )

    cuts = OrderedDict{Symbol, Vector{JuMP.ConstraintRef}}(
        :feasibility => JuMP.ConstraintRef[],
        :optimality => JuMP.ConstraintRef[],
    )

    log::Vector{String} = String[]
    timer::TimerOutput = TimerOutput()
end
DecomposedModel = DecomposedModel5

Base.show(io::IO, model::DecomposedModel) = print(io, "DecomposedModel [algorithm=Benders]: $(model.name)")

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

@kwdef struct DecomposedModel9 <: _DecompositionMainModule.AbstractDecomposedModel
    jump_model::JuMP.Model
    annotator::_DecompositionMainModule.AbstractGeneralAnnotator

    name::String = JuMP.name(jump_model)

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

    annotations = Dict{Symbol, Dict{Symbol, Any}}(
        :variables => Dict{Symbol, Any}(),
        :constraints => Dict{Symbol, Any}(),
    )

    # This allows for more than just `Benders.AbstractGeneralAttribute`, to allow (e.g.) `Solver.AbstractAttribute`.
    attributes::Vector{_DecompositionMainModule.DecompositionAttributeContainer} = _DecompositionMainModule.DecompositionAttributeContainer[]

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

    cuts = OrderedDict{Symbol, Vector{AbstractGeneralCut}}(
        :feasibility => AbstractGeneralCut[],
        :optimality => AbstractGeneralCut[],
    )

    log::Vector{String} = String[]
    timer::TimerOutput = TimerOutput()

    _cache::Dict{Symbol, Any} = Dict{Symbol, Any}()
    _defaults::Set{AbstractGeneralAttribute} = _ATTRIBUTE_DEFAULTS
end

# TODO: clean this up
DecomposedModel = DecomposedModel9

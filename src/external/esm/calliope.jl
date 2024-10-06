@kwdef struct Calliope <: AbstractExternalESM
    # TODO: these could keep version, specific necessary information, etc.
end

function _rows_with_entries(smcsc::SparseArrays.SparseMatrixCSC{Tv,Ti}, cols::AbstractVector{Ti}) where {Tv,Ti}
    nrows = size(smcsc, 1)
    row_has_one = falses(nrows)

    colptr = smcsc.colptr
    rowval = smcsc.rowval

    @inbounds @simd for j in cols
        idx_start = colptr[j]
        idx_end = colptr[j + 1] - 1
        row_has_one[rowval[idx_start:idx_end]] .= true
    end

    return findall(row_has_one)
end

function generate_annotation(model::Benders.DecomposedModel, ext_fw::Calliope)
    # TODO: make sure / check / warn on MILP models, since this, e.g., matches "available_flow_cap", which is not a design variable
    _split_vec(vec::Vector, L::Int) = [vec[round(Int, i*L) + 1:round(Int, (i+1)*L)] for i in 0:(length(vec) ÷ L - 1)]

    # Design variables are: flow_cap, link_flow_cap, source_cap, storage_cap, area_use
    #                       "cap" matches all cap-related of them
    #                       "a_use" matches area_use, without matching "source_use" (which is operational)
    REGEX_DESIGN = r".*(cap|a_use).*"

    # Get all variable names, and find the ones that match the design variables.
    var_names = JuMP.name.(model.lpmd.variables)
    vis_design = [i for i in axes(model.lpmd.A, 2) if !isnothing(match(REGEX_DESIGN, var_names[i]))]

    # Names of all constraints that link temporal "blocks".
    NAMES_TEMPORAL_LINKING_CONSTRAINTS = [
        "balance_supply_with_storage",  # i x T
        "balance_storage",              # i x T
        "ramping_up",                   # i x (T - 2)
        "ramping_down"                  # i x (T - 2)
    ]

    # Get all constraint names.
    con_names = JuMP.name.(model.lpmd.affine_constraints)

    # Prepare positive / negative entries of A.
    nzA = model.lpmd.A .!= 0
    posA = model.lpmd.A .> 0
    negA = model.lpmd.A .< 0

    # Prepare `single_var_con` for cache (will be converted to bounds).
    _srv = sort(model.lpmd.A.rowval)
    _single_var_con = Set{Int64}(
        _srv[i] for i in 2:(length(_srv) - 1)
        if (_srv[i] != _srv[i-1]) && (_srv[i] != _srv[i+1])
    )
    # Manually check the last constraint.
    if _srv[end] != _srv[end-1]
        push!(_single_var_con, _srv[end])
    end

    # Prepare information that is used in `model_from_lp`, to allow caching it.
    cache_model_from_lp = Dict(
        :nzA => nzA,
        :fnzc => nzA * collect(1:size(nzA, 2)),
        :single_var_con => _single_var_con,
    )

    # For each temporal linking constraint, find the corresponding variable indices, already split into blocks.
    vis_per_ntlc = Dict()
    for ntlc in NAMES_TEMPORAL_LINKING_CONSTRAINTS
        cis = [i for i in axes(model.lpmd.A, 1) if occursin(ntlc, con_names[i])]
        vis = findall(((sum(negA[cis, :]; dims=1) .== 1) .& (sum(posA[cis, :]; dims=1) .== 1))[1, :])

        if ntlc in ["balance_supply_with_storage", "balance_storage"]
            vis_per_ntlc[ntlc] = _split_vec(vis, model.T)
        elseif ntlc in ["ramping_up", "ramping_down"]
            vis_per_ntlc[ntlc] = vcat.(nothing, _split_vec(vis, model.T - 2), nothing)
        end
    end

    # TODO: the code above/below can be merged and made more efficient
    len_of_t_split = model.T ÷ model.nof_temporal_splits

    # Get all variables at the begin/end of each temporal block.
    vis_temporal = Int64[]
    for s in 1:model.nof_temporal_splits
        t_from = 1 + (s-1) * len_of_t_split
        t_to = s * len_of_t_split
    
        for (_, v) in vis_per_ntlc
            for el in v
                isnothing(el[t_from]) || push!(vis_temporal, el[t_from])
                isnothing(el[t_to]) || push!(vis_temporal, el[t_to])
            end
        end
    end

    # These variables will be present in main and should therefore NOT be included in the graph.
    set_vis_main = Set(vis_design)
    union!(set_vis_main, vis_temporal)
    vis_main = sort(collect(set_vis_main))

    # Get all constraints that do NOT contain variables that are not in the main-model.
    cis_main = findall((sum(nzA[:, [i for i in axes(model.lpmd.A, 2) if !(i in set_vis_main)]]; dims=2)[:, 1]) .== 0)

    # Construct the main-model.
    m_main = Benders.model_from_lp(model.lpmd, vis_main, cis_main; optimizer=model.f_opt_main(), cache=cache_model_from_lp)
    push!(model.models, m_main)
    push!(model.vis, vis_main)       
    push!(model.cis, cis_main)

    # Create a graph from the problem.
    adj_matrix = create_adjacency_matrix(model.lpmd.A, set_vis_main)
    g = Graphs.SimpleGraph(adj_matrix)
    cc = Graphs.connected_components(g)

    # We need row-wise access a lot of times later on, so transpose once, and reuse.
    # This is a big gain in performance. Copying is necessary, since `transpose` is lazy.
    nzAt = copy(nzA')

    # Create sub-model for each connected component.
    for component in cc
        if (length(component) > 1) || !(component[1] in set_vis_main)
            # These are the commands that we want to run, but substitute the transposed version AND use a specialized "find" function:
            #       `findall((sum(nzA[:, component]; dims=2) .!= 0)[:, 1])`
            #       `findall((sum(nzA[cis_in_component, :]; dims=1) .!= 0)[1, :])`
            # This results in around 100-500x speedup, depending on the size of the problem.
            cis_in_component = _rows_with_entries(nzA, component)
            vis_in_component = _rows_with_entries(nzAt, cis_in_component)
            
            m_sub = Benders.model_from_lp(model.lpmd, vis_in_component, cis_in_component; optimizer=model.f_opt_sub(), cache=cache_model_from_lp)
    
            push!(model.models, m_sub)
            push!(model.vis, vis_in_component)       
            push!(model.cis, cis_in_component)
        end
    end

    @info "Created decomposed models" n_sub = (length(model.models) - 1)

    return nothing
end

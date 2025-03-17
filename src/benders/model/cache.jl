_cache_build_nzA(model::Benders.DecomposedModel) = (model.lpmd.A .!= 0)
_cache_build_posA(model::Benders.DecomposedModel) = (model.lpmd.A .> 0)
_cache_build_negA(model::Benders.DecomposedModel) = (model.lpmd.A .< 0)
_cache_build_nzAt(model::Benders.DecomposedModel) = copy(cache_get(model, :nzA)')

function _cache_build_single_var_con(model::Benders.DecomposedModel)
    srv = sort(model.lpmd.A.rowval)
    single_var_con = Set{Int64}(srv[i] for i in 2:(length(srv)-1) if (srv[i] != srv[i-1]) && (srv[i] != srv[i+1]))

    # Manually check the last constraint.
    if srv[end] != srv[end-1]
        push!(single_var_con, srv[end])
    end

    return single_var_con
end

function _cache_build_fnzc(model::Benders.DecomposedModel)
    nzA = cache_get(model, :nzA)
    return nzA * collect(1:size(nzA, 2))
end

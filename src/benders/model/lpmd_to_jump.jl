function lpmd_to_jump(model::DecomposedModel, vis::Vector{Int64}, cis::Vector{Int64}; name::Union{Nothing, String} = nothing, optimizer)
    lpmd = model.lpmd

    cfg_direct_mode = get_attribute(model, Config.ModelDirectMode, :enable)
    cfg_debug = get_attribute(model, Config.ModelDebug, :enable)

    # Check if the model should be solved using a dual optimizer.
    dualize = false
    if startswith(name, "main")
        for attr in get_attributes(model, _ExtModSolver.DualizeModel)
            (attr.activate && attr.model == :main) || continue
            dualize = true
        end
    elseif startswith(name, "sub")
        for attr in get_attributes(model, _ExtModSolver.DualizeModel)
            (attr.activate && attr.model == :sub) || continue
            dualize = true
        end
    else
        @error "Unexpected model name" name
    end
    # TODO: we can use this name check above to stop requiring the "optimizer" argument, pulling it automatically

    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too

    jump_model = (
        if dualize
            # f_opt_dual = Dualization.dual_optimizer(() -> optimizer)
            # cfg_direct_mode ? JuMP.direct_model(f_opt_dual()) : JuMP.Model(f_opt_dual)

            # TODO: Is there a way to support dualization with direct models?
            if cfg_direct_mode
                @warn "Direct mode is not supported for dualized models, ignoring it" maxlog = 1
            end

            JuMP.Model(Dualization.dual_optimizer(() -> optimizer))
        else
            cfg_direct_mode ? JuMP.direct_model(optimizer) : JuMP.Model(() -> optimizer)
        end
    )
    
    JuMP.set_silent(jump_model)
    JuMP.set_string_names_on_creation(jump_model, cfg_debug)
    isnothing(name) || JuMP.set_attribute(jump_model, MOI.Name(), name)

    # Create variables, and set bounds.
    JuMP.@variable(jump_model, x[i = vis])
    for i in vis
        cfg_debug && JuMP.set_name(x[i], JuMP.name(lpmd.variables[i]))

        isfinite(lpmd.x_lower[i]) && JuMP.set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && JuMP.set_upper_bound(x[i], lpmd.x_upper[i])
    end

    # Get from cache:
    single_var_con = cache_get(model, :single_var_con)
    fnzc = cache_get(model, :fnzc)

    # Create constraints.
    add_ge = Vector{Int64}()
    add_lt = Vector{Int64}()
    add_eq = Vector{Int64}()
    for i in cis
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

    # Instead of:
    # ```julia
    # JuMP.@constraint(model, (@view lpmd.A[add_eq, vis]) * x.data .== lpmd.b_lower[add_eq])
    # JuMP.@constraint(model, (@view lpmd.A[add_ge, vis]) * x.data .>= lpmd.b_lower[add_ge])
    # JuMP.@constraint(model, (@view lpmd.A[add_lt, vis]) * x.data .<= lpmd.b_upper[add_lt])
    # ```
    #
    # we optimize to:
    # ```julia
    # x_data_t = collect(x.data')
    # A_t = copy(lpmd.A[:, vis]')
    # 
    # isempty(add_eq) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_eq], 1.0)' .== lpmd.b_lower[add_eq])
    # isempty(add_ge) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_ge], 1.0)' .>= lpmd.b_lower[add_ge])
    # isempty(add_lt) || JuMP.@constraint(model, mat_vec_scalar(x_data_t, A_t[:, add_lt], 1.0)' .<= lpmd.b_upper[add_lt])
    # ```
    # which again can be optimized by only transpose-copying once:

    x_data_t = collect(x.data')
    A_t = copy(lpmd.A[vcat(add_eq, add_ge, add_lt), vis]')  # TODO: this line takes > 50% of the total time
    
    idx_eq = 1:length(add_eq)
    idx_ge = length(add_eq)+1:length(add_eq)+length(add_ge)
    idx_lt = length(add_eq)+length(add_ge)+1:length(add_eq)+length(add_ge)+length(add_lt)   

    isempty(add_eq) || JuMP.@constraint(jump_model, mat_vec_scalar(x_data_t, A_t[:, idx_eq], 1.0)' .== lpmd.b_lower[add_eq])
    isempty(add_ge) || JuMP.@constraint(jump_model, mat_vec_scalar(x_data_t, A_t[:, idx_ge], 1.0)' .>= lpmd.b_lower[add_ge])
    isempty(add_lt) || JuMP.@constraint(jump_model, mat_vec_scalar(x_data_t, A_t[:, idx_lt], 1.0)' .<= lpmd.b_upper[add_lt])

    return jump_model
end

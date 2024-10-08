function lpmd_to_jump(model::DecomposedModel, vis::Vector{Int64}, cis::Vector{Int64}; name::Union{Nothing, String} = nothing, optimizer)
    lpmd = model.lpmd

    cfg_direct_mode = get_attribute(model, Config.ModelDirectMode, :enable)
    cfg_debug = get_attribute(model, Config.ModelDebug, :enable)

    # Check if the model should be solved using a dual optimizer.
    is_a_sub_model = false
    dualize = false
    if startswith(name, "main")
        for attr in get_attributes(model, _ExtModSolver.DualizeModel)
            (attr.activate && attr.model == :main) || continue
            dualize = true
        end
    elseif startswith(name, "sub")
        is_a_sub_model = true
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
            # if cfg_direct_mode
            #     @warn "Direct mode is not supported for dualized models, ignoring it" maxlog = 1
            # end

            # JuMP.Model(Dualization.dual_optimizer(() -> optimizer))
            JuMP.Model()
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

    # Set-up the linking/communication constraints to the main-problem.
    if is_a_sub_model
        shared_vis = sort(intersect(vis, model.vis[1]))
        sub_only_vis = setdiff(vis, shared_vis)
        
        JuMP.@variable(jump_model, y[i = shared_vis])
        JuMP.@constraint(jump_model, fix, x[shared_vis] .== y)

        jump_model[:obj] = JuMP.AffExpr(0.0)
        for vi in sub_only_vis
            JuMP.add_to_expression!(jump_model[:obj], lpmd.c[vi], x[vi])
        end
        
        # Construct objective / optimization sense.
        if has_attribute_type(model, Benders.CutTypeMISFSZ)
            # Modify the sub-problem into a feasibility only one, following the MIS-FSZ cut type.
            
            # The linking variable that get's "violated" from the main-model.
            JuMP.@variable(jump_model, θ_star)

            # Move objective into constraint.
            JuMP.@constraint(jump_model, limit_objective, jump_model[:obj] <= θ_star)
            JuMP.@objective(jump_model, Min, 0)
        else
            JuMP.@objective(jump_model, Min, jump_model[:obj])
        end
    end

    if dualize
        if cfg_direct_mode
            @warn "Direct mode is not supported for dualized models, ignoring it" maxlog = 1
        end

        dual_jump_model = JuMP.Model(() -> optimizer)

        JuMP.set_silent(dual_jump_model)
        JuMP.set_string_names_on_creation(dual_jump_model, cfg_debug)
        
        dual_problem = Dualization.DualProblem(JuMP.backend(dual_jump_model))
        Dualization.dualize(JuMP.backend(jump_model), dual_problem; variable_parameters = vcat(JuMP.index.(y.data), JuMP.index(θ_star)))

        # Mapping from initial variable indices to the dualized variables (parameters).
        dual_jump_model.ext[:dualization_var_to_vi] = Dict(
            JuMP.VariableRef(dual_jump_model, dual_problem.primal_dual_map.primal_parameter[JuMP.index(y[i])]) => i
            for i in shared_vis
        )

        # Pre-process the quadratic objective function.
        _params = keys(dual_jump_model.ext[:dualization_var_to_vi])
        e_q_obj = JuMP.objective_function(dual_jump_model)

        dual_jump_model.ext[:dualization_obj_base] = e_q_obj.aff
        dual_jump_model.ext[:dualization_obj_param] = [
            (term.first.a in _params) ? (term.first.a, term.first.b, term.second) : (term.first.b, term.first.a, term.second)
            for term in e_q_obj.terms
        ]

        # Handle MIS-FSZ cut type.
        if has_attribute_type(model, Benders.CutTypeMISFSZ)
            # Get the θ_star parameter in the dualized model.
            dual_jump_model.ext[:dualization_θ_star] = JuMP.VariableRef(dual_jump_model, dual_problem.primal_dual_map.primal_parameter[JuMP.index(θ_star)])

            # Use that to find the `π_0`.
            for term in JuMP.objective_function(dual_jump_model).terms
                if term.first.a == dual_jump_model.ext[:dualization_θ_star]
                    dual_jump_model.ext[:dualization_π_0] = term.first.b
                    break
                end
                if term.first.b == dual_jump_model.ext[:dualization_θ_star]
                    dual_jump_model.ext[:dualization_π_0] = term.first.a
                    break
                end
            end

            # TODO: instead of `dual_jump_model.ext[:dualization_π]`, just register it as `dual_jump_model[:π]`

            dual_jump_model.ext[:dualization_π] = [elem[2] for elem in dual_jump_model.ext[:dualization_obj_param] if elem[2] != dual_jump_model.ext[:dualization_θ_star]]

            # Construct helper variables to access the L1 norm of π.
            idx_π = eachindex(dual_jump_model.ext[:dualization_π])
            JuMP.@variable(dual_jump_model, π_pos[i = idx_π] >= 0)
            JuMP.@variable(dual_jump_model, π_neg[i = idx_π] >= 0)
            # JuMP.@constraint(dual_jump_model, cons_π_abs, π_pos .- π_neg .== dual_jump_model.ext[:dualization_π])
            JuMP.@constraint(dual_jump_model, cons_π_pos, π_pos .>= dual_jump_model.ext[:dualization_π])
            JuMP.@constraint(dual_jump_model, cons_π_neg, π_neg .>= -dual_jump_model.ext[:dualization_π])
            JuMP.@expression(dual_jump_model, l1_norm_π, sum(π_pos) + sum(π_neg))

            # Construct the normalization expression.
            # TODO NOTE WRITING: observe the "-1"/sign changing because the conic duality results in "<= 0" duals
            ω_0 = -1.0  # TODO: make this a parameter
            JuMP.@expression(dual_jump_model, expr_cglp_normalization, l1_norm_π + ω_0 * dual_jump_model.ext[:dualization_π_0])

            # Add the "CGLP normalization condition" (the "reduced" version as suggested).
            JuMP.@constraint(dual_jump_model, cons_cglp_normalization, expr_cglp_normalization == 1)
    
            # TODO WRITING NOTE: this uses conic duality, leading to a completely different model structure than the paper
            #                    this could have (bad?) implications for the cut generation, solving times, etc.
            #                    discuss / investigate this further (different dual formulations)
        end

        dual_jump_model.ext[:dualization_is_dualized] = true
        dual_jump_model.ext[:dualization_dual_problem] = dual_problem
        # dual_jump_model.ext[:dualization_primal_problem] = jump_model

        # TODO: remove this and account for it in cut calculation
        @assert unique(getindex.(dual_jump_model.ext[:dualization_obj_param], 3)) == [1.0]

        return dual_jump_model
    end

    return jump_model
end

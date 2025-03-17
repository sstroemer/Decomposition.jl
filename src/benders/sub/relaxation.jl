# TODO: use the "fix" constraint and inject the relaxation directly there

function apply!(model::Benders.DecomposedModel, attribute::Benders.Sub.RelaxationAll)
    @warn "`RelaxationAll` is currently untested ..."

    for i in 1:(length(model.models)-1)
        attribute.index == -1 || attribute.index == i || continue
        JuMP.relax_with_penalty!(Benders.sub(model; index = i); default = attribute.penalty)
    end

    # TODO: consider everything below ...

    # Get all variables that are part of the main-model, and find the constraints related to them.
    # idx_v = model.idx_model_vars[1]
    # idx_c_torelax = findall(sum(model.lpmd.A[:, idx_v] .!= 0; dims=2)[:, 1] .!= 0)

    # TODO: only relax the necessary constraints ...
    # relax_with_penalty!(bd_sub(model), Dict(linking_constraints .=> penalty))
    # merge!(attribute.penalty_map, relax_with_penalty!(m; default = attribute.penalty))

    return true
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Sub.RelaxationLinked)
    vis_main = model.vis[1]

    for i in 1:(length(model.models)-1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = Benders.sub(model; index = i)
        vis_main_in_sub = [vi for vi in vis_main if vi in m_sub[:y].axes[1]]

        JuMP.@variable(m_sub, z_pos[i = vis_main_in_sub], lower_bound = 0)
        JuMP.@variable(m_sub, z_neg[i = vis_main_in_sub], lower_bound = 0)

        all_con = JuMP.all_constraints(m_sub; include_variable_in_set_constraints = false)

        exp_penalty = JuMP.AffExpr(0.0)
        for vi in vis_main_in_sub
            coeffs = JuMP.normalized_coefficient.(all_con, m_sub[:y][vi])
            JuMP.set_normalized_coefficient.(all_con, z_pos[vi], coeffs)
            JuMP.set_normalized_coefficient.(all_con, z_neg[vi], -coeffs)

            JuMP.add_to_expression!(exp_penalty, z_pos[vi], attribute.penalty)
            JuMP.add_to_expression!(exp_penalty, z_neg[vi], attribute.penalty)
        end

        JuMP.@objective(m_sub, Min, JuMP.objective_function(m_sub) + exp_penalty)
    end

    return true
end

function apply!(model::Benders.DecomposedModel, attribute::Benders.Sub.RelaxationRegex)
    vn_main = JuMP.name.(model.lpmd.variables[model.vis[1]])
    vis_filtered_main = [vi for (vi, vn) in zip(model.vis[1], vn_main) if !isnothing(match(attribute.regex, vn))]

    for i in 1:(length(model.models)-1)
        attribute.index == -1 || attribute.index == i || continue

        m_sub = sub(model; index = i)
        vis_main_in_sub = [vi for vi in vis_filtered_main if vi in m_sub[:y].axes[1]]

        JuMP.@variable(m_sub, z_pos[i = vis_main_in_sub], lower_bound = 0)
        JuMP.@variable(m_sub, z_neg[i = vis_main_in_sub], lower_bound = 0)

        all_con = JuMP.all_constraints(m_sub; include_variable_in_set_constraints = false)

        exp_penalty = JuMP.AffExpr(0.0)
        for vi in vis_main_in_sub
            coeffs = JuMP.normalized_coefficient.(all_con, m_sub[:y][vi])
            JuMP.set_normalized_coefficient.(all_con, z_pos[vi], coeffs)
            JuMP.set_normalized_coefficient.(all_con, z_neg[vi], -coeffs)

            JuMP.add_to_expression!(exp_penalty, z_pos[vi], attribute.penalty)
            JuMP.add_to_expression!(exp_penalty, z_neg[vi], attribute.penalty)
        end

        JuMP.@objective(m_sub, Min, JuMP.objective_function(m_sub) + exp_penalty)
    end

    return true
end

function execute!(model::Benders.DecomposedModel, query::Benders.Query.Relaxation)
    # TODO: account for different previous ways to ensure feasibility (by checking the list of attributes/modifications)
    violation = sum(JuMP.value.(sub(model)[:z_pos])) + sum(JuMP.value.(sub(model)[:z_neg]))
    return isapprox(violation, 0.0; atol = 1e-6)
end

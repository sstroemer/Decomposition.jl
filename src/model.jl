abstract type DecompositionAttribute end
abstract type DecompositionQuery end

@kwdef struct DecomposedModel4
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

    stats = Dict{Symbol, Any}(
        :iteration => 0,
        :lower_bound => -Inf,
        :upper_bound => Inf,
        :gap_rel => Inf,
        :gap_abs => Inf,
    )     # TODO: transform into a struct, that tracks stats for each iteration (using a "inc_iter" function)

    attributes::Vector{DecompositionAttribute} = DecompositionAttribute[]

    cuts = Dict{Symbol, Vector{ConstraintRef}}(
        :feasibility => ConstraintRef[],
        :optimality => ConstraintRef[],
    )  # TODO: track which cut is created by which iteration (inside the stats struct)
end
DecomposedModel = DecomposedModel4

function modify(::DecomposedModel, ::DecompositionAttribute)
    @error "Not implemented"
end

attach_solver(model::JuMP.Model, solver) = set_optimizer(model, solver)  # TODO: use this to attach "BD" (or others)
solve!(model::JuMP.Model) = optimize!(model)

function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64})
    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too
    # TODO: transform single variable constraints into bounds

    model = direct_model(Gurobi.Optimizer())

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
    n = size(A, 2)
    ama = SparseArrays.spzeros(Int, n, n)

    # @infiltrate
    # nzA = A .!= 0
    # rvs = rowvals(nzA)
    # for vi in 1:size(A, 2)
    #     rows = rvs[nzrange(nzA, vi)]

    #     for row in rows
    #         for var in A[row, :].nzind
    #             var == vi && continue
    #             ama[vi, var] = 1
    #         end
    #     end
    # end

    # # Variable i is contained in
    # rows = rowvals(A)[nzrange(nzA, vi)]
    # # which contains the following total variables
    # vars = A[853, :].nzind

    for i in 1:size(A, 1)
        nodes = A[i, :].nzind
        for j in 1:length(nodes)
            u = nodes[j]
            nodes[j] in rm && continue
            for k in j+1:length(nodes)
                v = nodes[k]
                nodes[k] in rm && continue
                ama[u, v] = 1
                # ama[v, u] = 1
            end
        end
    end

    return LinearAlgebra.Symmetric(ama)
end
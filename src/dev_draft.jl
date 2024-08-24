using JuMP
import HiGHS
import Statistics: mean
using Chairmarks

import Printf
function _print_iteration(k, args...)
    f(x) = Printf.@sprintf("%11.3e", x)
    println(lpad(k, 8), " │ ", join(f.(args), " │ "))
    return
end


# SCALING ================================================================

# mtype = "normal"
# model = mtype == "normal" ? Model(HiGHS.Optimizer) : direct_model(HiGHS.Optimizer())

model = read_from_file("ns070.MPS"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
set_optimizer(model, HiGHS.Optimizer)   # 2.7704847864e+05

mmd = lp_matrix_data(model)
A = mmd.A

println("       k |      min(c) |      avg(c) |      max(c) |    range(c) |             |      min(A) |      avg(A) |      max(A) |    range(A)")

range_c = extrema(abs.(filter(!=(0), mmd.c)))
range_A = extrema(abs.(mmd.A.nzval))
c_avg = mean(abs.(filter(!=(0), mmd.c)))
A_avg = mean(abs.(mmd.A.nzval))
_print_iteration("init", range_c[1], c_avg, range_c[2], range_c[2]/range_c[1], Inf, range_A[1], A_avg, range_A[2], range_A[2]/range_A[1])

# TODO: some of the following steps require scaling `b` or `lb`/`ub` as well - fix this!
mmd.A.nzval[:] /= range_A[1]
mmd.c ./= range_c[1]

range_c = extrema(abs.(filter(!=(0), mmd.c)))
range_A = extrema(abs.(mmd.A.nzval))
c_avg = mean(abs.(filter(!=(0), mmd.c)))
A_avg = mean(abs.(mmd.A.nzval))
_print_iteration("basic", range_c[1], c_avg, range_c[2], range_c[2]/range_c[1], Inf, range_A[1], A_avg, range_A[2], range_A[2]/range_A[1])

for k in 1:10
    range_c = extrema(abs.(filter(!=(0), mmd.c)))
    range_A = extrema(abs.(mmd.A.nzval))

    c_avg = mean(abs.(filter(!=(0), mmd.c)))
    A_avg = mean(abs.(mmd.A.nzval))

    _print_iteration(k, range_c[1], c_avg, range_c[2], range_c[2]/range_c[1], Inf, range_A[1], A_avg, range_A[2], range_A[2]/range_A[1])

    for j in 1:size(mmd.A, 2)
        c_j = abs(mmd.c[j])
        A_view = @view mmd.A.nzval[mmd.A.colptr[j]:(mmd.A.colptr[j + 1] - 1)]
        A_j_range = extrema(abs.(A_view); init = (+Inf, -Inf))

        A_j_range[1] == +Inf && continue    # TODO: prevent this and just scale `c` instead

        min_scale_c = 1.0 / (range_c[2] / c_j)
        max_scale_c = 1.0 / (range_c[1] / c_j)

        min_scale_A = 1.0 / (range_A[2] / maximum(A_j_range))
        max_scale_A = 1.0 / (range_A[1] / minimum(A_j_range))

        # NOTE: we could scale the "extreme" rows of this col, to influence scaling
        
        if (c_j / range_c[1]) > (range_c[2] / range_c[1]) * 0.25
            # high c => x_new ~ 1000 * x
        end

        mmd.c[j] /= scale
        mmd.A[:, j] ./= scale   # TODO: do that using the nzval accessor
    end
end


# =====================================================================
# =====================================================================
# =====================================================================
# =====================================================================


@kwdef struct DecompDataDict
    _internal::Dict = Dict{Any, Any}(
        :annotations => Dict{Any, Any}()
    )
end
Base.setproperty!(x::DecompDataDict, s::Symbol, v) = getfield(x, :_internal)[s] = v
Base.getproperty(x::DecompDataDict, s::Symbol) = getfield(x, :_internal)[s]
DecompositionData = DecompDataDict
_djl(m::Model) = m.ext[:decomposition]

abstract type DecompositionAttribute end
# TODO: Keep a list of all "attributes" (which are just sub-types that describe modifications, and may contain more data (e.g., the sub-model index they were applied to))

@kwdef struct DecomposedModel12
    monolithic::JuMP.Model

    lpmd::JuMP.LPMatrixData

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

    cuts = Vector{Any}()  # TODO: track which cut is created by which iteration (inside the stats struct)
end
DecomposedModel = DecomposedModel12

bd_main(model::DecomposedModel) = model.models[1]
bd_sub(model::DecomposedModel; index::Int = 1) = model.models[1 + index]
bd_subs(model::DecomposedModel) = model.models[2:end]

function add_annotation!(model::DecomposedModel, v::VariableRef, annotation::Symbol)
    if !haskey(model.annotations[:variables], v)
        model.annotations[:variables][v] = Set{Symbol}()
    end
    push!(model.annotations[:variables][v], annotation)

    return nothing
end

function has_annotation(model::DecomposedModel, v::VariableRef, annotation::Symbol)
    haskey(model.annotations[:variables], v) || return false
    return annotation in model.annotations[:variables][v]
end

function bd_derive_constraint_annotations(model::DecomposedModel)
    # TODO: do not derive annotations for constraints that were manually annotated

end

function bd_decompose(model::DecomposedModel)
    if !isempty(model.models)
        @warn "Overwriting an existing decomposition"
        empty!(model.models)
        empty!(model.reference_maps)
    end

    A = model.monolithic.ext[:lp_matrix_data].A
    v_main = BitVector(has_annotation(model, v, :design) for v in all_variables(model.monolithic))
    v_sub = .~v_main
    nzA = (A .!= 0)
    has_main_vars = nzA * v_main .!= 0
    has_sub_vars = nzA * v_sub .!= 0

    f_filter_main(c) = has_main_vars[index(c).value] && !has_sub_vars[index(c).value]
    f_filter_sub(c) = has_sub_vars[index(c).value]

    m_main, rm_main = copy_model(model.monolithic; filter_constraints=f_filter_main)
    push!(model.models, m_main)
    push!(model.reference_maps, rm_main)

    # TODO: Support more than 1 sub problem
    m_sub, rm_sub = copy_model(model.monolithic)
    push!(model.models, m_sub)
    push!(model.reference_maps, rm_sub)

    # Change interval constraints to bounds, since they might interfere with fixing variables.
    for c_index in all_constraints(m_sub, VariableRef, MOI.Interval{Float64})
        constraint = constraint_object(c_index)
        delete(m_sub, c_index)
        set_lower_bound(constraint.func, constraint.set.lower)
        set_upper_bound(constraint.func, constraint.set.upper)
    end

    delvars = VariableRef[]
    for v in all_variables(model.monolithic)
        if v_main[index(v).value]
            # This is a main variable.
            model.decomposition_maps[rm_main[v]] = [rm_sub[v]]  # TODO: account for more than 1 sub
            set_objective_coefficient(m_sub, rm_sub[v], 0)

            # Delete all constraints that are only attached to this variable.
            # has_lower_bound(rm_sub[v]) && delete_lower_bound(rm_sub[v])
            # has_upper_bound(rm_sub[v]) && delete_upper_bound(rm_sub[v])
            # `VariableIndex`-in-`MathOptInterface.Interval{Float64}
            continue
        else
            # This is NOT a main variable.
            push!(delvars, rm_main[v])
            # delete(m_main, rm_main[v])
            # TODO: delete!(rm_main.index_map.var_map, v)
        end
    end
    delete(m_main, delvars)

    # Create sub-model estimation variable(s).
    # TODO: allow `θ[s = 1:S]`
    @variable(m_main, θ >= 0)

    # Create objective function.
    @expression(m_main, exp_objective_original, objective_function(m_main))
    @expression(m_main, exp_objective_θ, m_main[:θ])
    @objective(m_main, Min, exp_objective_θ + exp_objective_original)

    return nothing
end

function bd_modify_main_ensure_var_bounds(model::DecomposedModel; bound::Float64 = 1e8)
    # TODO: make sure this DOES NOT happen for θ -> solvable if "all_variables" is kept in struct
    for v in all_variables(bd_main(model))
        has_lower_bound(v) || set_lower_bound(v, -bound)
        has_upper_bound(v) || set_upper_bound(v, bound)
    end

    return nothing
end

function bd_modify_sub_ensure_feasibility(model::DecomposedModel; index::Int64 = -1, penalty::Float64 = 1e7)
    # TODO: handle more than 1 sub-model, and use the index (-1 = all, 0 = error, 1 = first, ...)

    linking_constraints = Vector{ConstraintRef}()   # TODO: this should be stored inside the deomcposedmodel
    for var in all_variables(bd_main(model))
        haskey(model.decomposition_maps, var) || continue   # TODO: remove this, as soon as θ is not part of the decomposition map anymore
        sub_var = model.decomposition_maps[var][1]::JuMP.VariableRef  # TODO: handle index of sub-models
        for con in all_constraints(bd_sub(model; index=1); include_variable_in_set_constraints = false)
            # TODO: check `include_variable_in_set_constraints` needed?
            # TODO: read performance note in the docs, and use a function barrier
            if normalized_coefficient(con, sub_var) != 0
                push!(linking_constraints, con)
            end
        end
    end

    penalty_map = relax_with_penalty!(bd_sub(model), Dict(linking_constraints .=> penalty))
    
    return nothing
end

attach_solver(model::JuMP.Model, solver) = set_optimizer(model, solver)  # TODO: use this to attach "BD" (or others)

solve!(model::JuMP.Model) = optimize!(model)

function bd_insert_cuts!(model::DecomposedModel)
    # TODO: account for more than 1 sub-model

    exp_cut = AffExpr(objective_value(bd_sub(model; index=1)))
    for var in all_variables(bd_main(model))
        haskey(model.decomposition_maps, var) || continue
        π = reduced_cost(model.decomposition_maps[var][1])
        add_to_expression!(exp_cut, π, var - value(var))
    end

    push!(model.cuts, @constraint(bd_main(model), bd_main(model)[:θ] >= exp_cut))

    return nothing
end

function iterate!(model::DecomposedModel)
    # TODO: make that a "setting"
    set_silent(bd_main(model))
    set_silent.(bd_subs(model))

    # TODO: account for more than 1 sub-model

    if model.stats[:iteration] == 1
        has_lower_bound(bd_main(model)[:θ]) && delete_lower_bound(bd_main(model)[:θ])
        has_upper_bound(bd_main(model)[:θ]) && delete_upper_bound(bd_main(model)[:θ])
    end

    # Solve main-model.
    solve!(bd_main(model))

    # Update lower bound.
    model.stats[:lower_bound] = objective_value(bd_main(model))

    # Update sub-models based on new values from main-model.
    for var in all_variables(bd_main(model))
        haskey(model.decomposition_maps, var) || continue
        for var_sub in model.decomposition_maps[var]
            fix(var_sub, value(var); force=true)
        end
    end

    # Solve sub-models.
    for m in bd_subs(model)
        solve!(m)
    end

    # Update upper bound.
    ub_candidate = (
        value(bd_main(model)[:exp_objective_original]) +
        sum(objective_value(m) for m in bd_subs(model))
    )
    if ub_candidate < model.stats[:upper_bound]
        model.stats[:upper_bound] = ub_candidate
        # TODO: update best solution
    end

    # Update stats.
    # TODO: Move this into a function that is bound to the "stats" struct
    model.stats[:gap_abs] = model.stats[:upper_bound] - model.stats[:lower_bound]
    model.stats[:gap_rel] = model.stats[:gap_abs] / abs(model.stats[:upper_bound])  # TODO: check, use ub or lb here? (gurobi)
    
    # Update cuts.
    bd_insert_cuts!(model)

    # TODO: Move the printing into the next_iter function
    _print_iteration(model.stats[:iteration], model.stats[:lower_bound], model.stats[:upper_bound], model.stats[:gap_rel])
    model.stats[:iteration] += 1

    return nothing
end


orig_model = read_from_file("uc_070_1y.mps"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
orig_model = read_from_file("ns070.MPS"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
_om, _ = copy_model(orig_model)

# TODO: check that variables and constraints are 1:N indexed, since we currently blindly trust this (e.g., for efficient matrix ops)
# TODO: check (or later handle) that the model is a minimization problem

model = DecomposedModel(; monolithic = orig_model, lpmd = lp_matrix_data(orig_model))

design_vs = ["flow_cap", "link_flow_cap", "source_cap", "storage_cap", "area_use"]
vs = all_variables(orig_model)
vv = [v for v in vs if any(occursin(el, name(v)) for el in design_vs)]

for var in vv
    add_annotation!(model, var, :design)
end

# TODO: transform `bd_decompose` into `Benders.decompose!(...)` with a proper submodule `Benders`
bd_decompose(model)
bd_modify_main_ensure_var_bounds(model; bound = 1e7)

attach_solver(bd_main(model), HiGHS.Optimizer)
set_attribute(bd_main(model), "solver", "ipm")
set_attribute(bd_main(model), "run_crossover", "off")
attach_solver(bd_sub(model; index=1), HiGHS.Optimizer)

bd_modify_sub_ensure_feasibility(model; penalty=1e5)

iterate!(model)


# NOTE FOR TOMORROW:
# are the "interval" constraints copied correctly? if yes, then the "set bound" should fail, right?




# TODO: given a linked constraint x + 3y_1 - 2y_2 <= b, we can determine "bounds" for x, by looking at lower bounds of y_1 and upper bounds of y_2

# TODO: check sense, integer, binary, for Min, none, none


lpmd = model.monolithic.ext[:lp_matrix_data]




function model_from_lp(lpmd::JuMP.LPMatrixData, idx_v::Vector{Int64}, idx_c::Vector{Int64})
    # TODO: allow `::Base.OneTo{Int64}` instead of `::Vector{Int64}` too
    # TODO: transform single variable constraints into bounds

    model = direct_model(HiGHS.Optimizer())

    # Create variables, and set bounds.
    @variable(model, x[i = idx_v])
    for i in eachindex(x)
        isfinite(lpmd.x_lower[i]) && set_lower_bound(x[i], lpmd.x_lower[i])
        isfinite(lpmd.x_upper[i]) && set_upper_bound(x[i], lpmd.x_upper[i])
    end

    # Create constraints.
    for i in idx_c
        if lpmd.b_lower[i] == lpmd.b_upper[i]
            # If `lb == ub`, then we know both have to be finite.
            @constraint(model, lpmd.A[i, idx_v]' * x.data == lpmd.b_lower[i])
            continue
        end

        isfinite(lpmd.b_lower[i]) && @constraint(model, lpmd.A[i, idx_v]' * x.data >= lpmd.b_lower[i])
        isfinite(lpmd.b_upper[i]) && @constraint(model, lpmd.A[i, idx_v]' * x.data <= lpmd.b_upper[i])
    end

    return model
end

"""

main: obj, con
sub: pure, full (includes main costs)
"""
struct BD_MainObjectiveObj <: DecompositionAttribute end
struct BD_MainObjectiveCon <: DecompositionAttribute end

struct BD_SubObjectivePure <: DecompositionAttribute; index::Int64; end
struct BD_SubObjectiveFull <: DecompositionAttribute; index::Int64; end

function bd_modify(model::DecomposedModel, attribute::DecompositionAttribute)
   @error "Not implemented"
end

function bd_modify(model::DecomposedModel, ::BD_MainObjectiveObj)
    m = bd_main(model)
    idx_v = model.idx_model_vars[1]
    
    @objective(m, Min, model.lpmd.c[idx_v]' * m[:x].data + model.lpmd.c_offset)

    return nothing
end

function bd_modify(model::DecomposedModel, attribute::BD_SubObjectivePure)
    m = bd_sub(model; index=attribute.index)
    idx_v = model.idx_model_vars[1 + attribute.index]

    obj = AffExpr(0.0)
    for i in idx_v
        (i in model.idx_model_vars[1]) && continue
        add_to_expression!(obj, model.lpmd.c[i], m[:x][i])
    end
    @objective(m, Min, obj)

    return nothing
end




using JuMP
import HiGHS
import Statistics: mean
using Chairmarks

import Printf
function _print_iteration(k, args...)
    f(x) = Printf.@sprintf("%11.3e", x)
    println(lpad(k, 8), " │ ", join(f.(args), " │ "))
    return
end

include("draft/model.jl")
include("draft/benders.jl")

# --------------------------------------------------------------------

orig_model = read_from_file("ns070.MPS"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
model = DecomposedModel(; monolithic = orig_model, lpmd = lp_matrix_data(orig_model))

design_vs = ["flow_cap", "link_flow_cap", "source_cap", "storage_cap", "area_use"]
design_vs_idx = [index(v).value for v in all_variables(orig_model) if any(occursin(el, name(v)) for el in design_vs)]

v_main = BitVector(any(occursin(el, name(v)) for el in design_vs) for v in all_variables(orig_model))
v_sub = .~v_main
nzA = (model.lpmd.A .!= 0)
has_main_vars = nzA * v_main .!= 0
has_sub_vars = nzA * v_sub .!= 0

idx_main_vars = design_vs_idx
idx_main_cons = findall(has_main_vars .& (.~has_sub_vars))

# --------------------------------------------------------------------

m_main = model_from_lp(model.lpmd, idx_main_vars, idx_main_cons)
m_sub = model_from_lp(model.lpmd, collect(axes(model.lpmd.A, 2)), collect(i for i in axes(model.lpmd.A, 1) if i ∉ idx_main_cons))

push!(model.models, m_main)
push!(model.models, m_sub)

push!(model.idx_model_vars, idx_main_vars)
push!(model.idx_model_vars, collect(axes(model.lpmd.A, 2)))

push!(model.idx_model_cons, idx_main_cons)
push!(model.idx_model_cons, collect(i for i in axes(model.lpmd.A, 1) if i ∉ idx_main_cons))

# --------------------------------------------------------------------

bd_modify(model, BD_MainObjectiveObj())
bd_modify(model, BD_SubObjectivePure(1))

bd_modify(model, BD_SubEnsureFeasibility(1, 1e6))

optimize!(bd_main(model))

sol = value.(bd_main(model)[:x])
for s in 1:(length(model.models) - 1), i in eachindex(sol)
    m = bd_sub(model; index=s)
    fix(m[:x][i], sol[i]; force=true)
end

optimize!(bd_sub(model; index=1))


exp_cut = AffExpr(objective_value(bd_sub(model; index=1)))
for i in eachindex(sol)
    π = reduced_cost(bd_sub(model; index=s)[:x][i])
    add_to_expression!(exp_cut, π, bd_main(model)[:x][i] - value(var))
end

push!(model.cuts, @constraint(bd_main(model), bd_main(model)[:θ] >= exp_cut))



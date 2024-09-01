# TODO: check that variables and constraints are 1:N indexed, since we currently blindly trust this (e.g., for efficient matrix ops)
# TODO: check (or later handle) that the model is a minimization problem
# TODO: given a linked constraint x + 3y_1 - 2y_2 <= b, we can determine "bounds" for x, by looking at lower bounds of y_1 and upper bounds of y_2
# TODO: check sense, integer, binary, for Min, none, none


# import JuMP
import HiGHS, Gurobi
using Decomposition

# const GRB_ENV = Gurobi.Env()

# using Logging
# SWITCH_TO_DEBUG = false
# if SWITCH_TO_DEBUG
#     global_logger(ConsoleLogger(stderr, Logging.Debug))
# else
#     global_logger(ConsoleLogger(stderr, Logging.Info))
# end

T = 120
nof_temporal_splits = 5
jump_model = jump_model_from_file("national_scale_$T.mps")

# Decomposition.JuMP.set_optimizer(jump_model, Gurobi.Optimizer)
# Decomposition.JuMP.optimize!(jump_model)

model = Benders.DecomposedModel(; jump_model, T, nof_temporal_splits, f_opt = Gurobi.Optimizer)
generate_annotation(model, Calliope())

modify(model, Solver.AlgorithmIPM(model = :main))

modify(model, Solver.ExtractDualRay(model = :sub))
modify(model, Benders.Main.FeasibilityCutTypeMulti())
modify(model, Benders.Main.OptimalityCutTypeMulti())

modify(model, Benders.Main.ObjectiveDefault())
modify(model, Benders.Sub.ObjectiveSelf())

modify(model, Benders.Main.VirtualSoftBounds(0.0, 1e6))

# modify(model, Benders.Sub.RelaxationLinked(-1, 1e6))

modify(model, Benders.Termination.Stop(opt_gap_rel = 1e-2))

for k in 1:100
    iterate!(model) && break
end

save(model)
# TODO: bd_query(model, BD_SubFeasibility())







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++
# TODO: HIGH PRIO
# READ::::: https://discourse.julialang.org/t/simplifying-lp-for-repeated-resolving/97528/10
# ==> warmstarting sub problems ==> Lene

# Q: Experimented with switching to single cut later on? (seems to be "smoother" in converging after being close)
# Q: Am I dumb? (discourse)
# Q: "Approx" pcores - how?
# Q: subproblems only with barrier?

# Cluster by complexity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~++




include("visualization.jl")


_legend_names = ["feas. cuts: on", "feas. cuts: off", "feas. cuts: storage"]
compare(
    ["a064ac2", "e757431", "4cd40b2"],
    [
        trace_rel_gap,
        (mi) -> (
            x=0:length(mi["history"]) - 1,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Some title goes here",
        yaxis_type="log",
        xaxis_title="Iteration",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="01_rel_gap_iter"
)

_legend_names = ["feas. cuts: on", "feas. cuts: off", "feas. cuts: storage"]
compare(
    ["a064ac2", "e757431", "4cd40b2"],
    [
        (mi) -> trace_bounds(mi)[1],
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Some title goes here",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Lower bound",
        xaxis_range=[0, 30],
        yaxis_range=log10.([1e3, 3e5]),
    );
    html = "nothing"
)























# --------------------------------------------------------------------
# REL GAP over ITERATION
# --------------------------------------------------------------------
compare(
    ["88725ff", "a311d66", "129ea42", "ca16bda", "cb707c9"],
    [
        trace_rel_gap,
        (mi) -> (
            x=0:length(mi["history"]) - 1,
            name="$(mi["inputs"]["splits"]) splits",
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Iteration",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="01_rel_gap_iter"
)

# --------------------------------------------------------------------
# REL GAP over TIME
# --------------------------------------------------------------------
compare(
    ["88725ff", "a311d66", "129ea42", "ca16bda", "cb707c9"],
    [
        trace_rel_gap,
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name="$(mi["inputs"]["splits"]) splits",
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="02_rel_gap_time"
)

# --------------------------------------------------------------------
# LOWER BOUND over TIME
# --------------------------------------------------------------------
compare(
    ["88725ff", "a311d66", "129ea42", "ca16bda", "cb707c9"],
    [
        (mi) -> trace_bounds(mi)[1],
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name="$(mi["inputs"]["splits"]) splits",
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Lower bound",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Lower bound",
        xaxis_range=[0, 20],
        yaxis_range=log10.([2e5, 13e6]),
    );
    additional_traces=[
        scatter(x = [0, 60], y = repeat([1.1148105940e+07], 2), mode="lines", name="true solution", line=attr(color="black", width=3, dash="dot")),
        scatter(x = [3.08, 3.08], y = [0, 15e6], mode="lines", name="from gurobi", line=attr(color="black", width=3, dash="dot")),
    ],
    html="03_lb_time"
)

# --------------------------------------------------------------------
# (T:8760) LOWER BOUND over TIME
# --------------------------------------------------------------------
compare(
    ["5b54bec"],
    [
        (mi) -> trace_bounds(mi)[1],
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name="$(mi["inputs"]["splits"]) splits",
        )
    ],
    Layout(
        title="Calliope - National Scale (1Y) - Lower bound",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Lower bound",
        xaxis_range=[0, 25],
        yaxis_range=log10.([1e5, 25e6]),
    );
    additional_traces=[
        scatter(x = [0, 60], y = repeat([2.238932354e+07], 2), mode="lines", name="true solution", line=attr(color="black", width=3, dash="dot")),
        scatter(x = [18.43, 18.43], y = [0, 25e6], mode="lines", name="from gurobi", line=attr(color="black", width=3, dash="dot")),
    ],
    html="04_lb_time_full_year"
)

# --------------------------------------------------------------------
# (select:split:48) REL GAP over ITERATION
# --------------------------------------------------------------------
_legend_names = ["feasibility (on demand)", "feasibility (always)", "relax (full)", "relax (storage)"]
compare(
    ["129ea42", "9b0ff83", "0a91362", "9c10568"],
    [
        trace_rel_gap,
        (mi) -> (
            x=0:length(mi["history"]) - 1,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Iteration",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="05_rel_gap_iter_feasibility"
)

# --------------------------------------------------------------------
# (select:split:48) REL GAP over TIME
# --------------------------------------------------------------------
_legend_names = ["feasibility (on demand)", "feasibility (always)", "relax (full)", "relax (storage)"]
compare(
    ["129ea42", "9b0ff83", "0a91362", "9c10568"],
    [
        trace_rel_gap,
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="06_rel_gap_time_feasibility"
)

# --------------------------------------------------------------------
# (select:split:48:feasibility:ondemand) REL GAP over ITERATION
# --------------------------------------------------------------------
_legend_names = ["default (dual simplex)", "sub: predual", "sub: predual + numericfocus=3", "sub: primal simplex", "sub: lpwarmstart=2", "sub: lpwarmstart=0", "sub: lpwarmstart cond."]
compare(
    ["129ea42", "a0cb442", "ff144d9", "826d4ec", "4401e47", "e497a03", "88a8421"],
    [
        trace_rel_gap,
        (mi) -> (
            x=0:length(mi["history"]) - 1,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Iteration",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="07_rel_gap_iter_solveroptions"
)

# --------------------------------------------------------------------
# (select:split:48:feasibility:ondemand) REL GAP over TIME
# --------------------------------------------------------------------
_legend_names = ["default (dual simplex)", "sub: predual", "sub: predual + numericfocus=3", "sub: primal simplex", "sub: lpwarmstart=2", "sub: lpwarmstart=0", "sub: lpwarmstart cond."]
compare(
    ["129ea42", "a0cb442", "ff144d9", "826d4ec", "4401e47", "e497a03", "88a8421"],
    [
        trace_rel_gap,
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Calliope - National Scale (0.5Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="08_rel_gap_iter_solveroptions"
)

# --------------------------------------------------------------------
# (T:2184:inner:highs) REL GAP over TIME
# --------------------------------------------------------------------
_legend_names = ["gurobi + gurobi", "gurobi + highs"]
compare(
    ["20624ac", "7a47a37"],
    [
        trace_rel_gap,
        (mi) -> (
            x=cumsum(it["time"]["wall"] for it in mi["history"]) ./ 1e9,
            name=popfirst!(_legend_names),
        )
    ],
    Layout(
        title="Calliope - National Scale (0.25Y) - Relative gap",
        yaxis_type="log",
        xaxis_title="Time [s]",
        yaxis_title="Relative gap",
        yaxis_range=log10.([1e-12, 5.0]),
    );
    html="09_rel_gap_iter_solvers"
)



different algorithms, different settings, tested with different solvers
=> so: what works with OS solvers, what do they need, etc.

handbook / guide / clustered overview
step by step overview of techniques, parameters, algorithms, etc.

"toy" examples for each "introduced change"


document the process




prefixes = ["asd", "5sd", "asdf"]

_(549bc21|549bc23|549b11)\.djl\.json$



readdir("out")

model_info = JSON3.read("out/national_scale_2184_549bc21.djl.json", Dict; allow_inf=true)

max_iter = length(model_info["history"]) - 1

lb = Vector{Float64}([model_info["history"][1]["lower_bound"]])
ub = Vector{Float64}([model_info["history"][1]["upper_bound"]])
for entry in model_info["history"][2:end]
    push!(lb, max(lb[end], entry["lower_bound"]))
    push!(ub, min(ub[end], entry["upper_bound"]))
end

plot([
    scatter(x=0:max_iter, y=rel_gap.(lb, ub), mode="lines", line_shape="hv")
])





_m = Model(HiGHS.Optimizer)

@variable(_m, x[1:2] >= 0)
@variable(_m, y[1:2])

fix.(y, 1.0; force=true)

c1 = @constraint(_m, c1, x[1] <= y[1])
c2 = @constraint(_m, c2, sum(x) >= 3)
c3 = @constraint(_m, c3, x[2] <= y[2])

@objective(_m, Min, -sum(x))

# for gurobi:
set_attribute(_m, "InfUnbdInfo", 1)
set_attribute(_m, "DualReductions", 0)
# for highs:
set_attribute(_m, "presolve", "off")

optimize!(_m)

primal_status(_m) == MOI.NO_SOLUTION
dual_status(_m) == MOI.NO_SOLUTION # => need presolve off

dual_status(_m) == MOI.INFEASIBILITY_CERTIFICATE    # => feasibility cut

ζ = dual.(FixRef.(y))
b = fix_value.(y)


@constraint(_m, ζ' * (y .- b) <= 0)
@constraint(_m, ζ' * y <= 0)


dual.([c1, c2, c3])


c = all_constraints(_m; include_variable_in_set_constraints = true)

dual.(c)




unset_silent(bd_main(model))
mwe = bd_main(model)
write_to_file(mwe, "mwe.mps.bz2"; format = JuMP.MOI.FileFormats.FORMAT_MPS)

mwe = read_from_file("mwe.mps"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
set_optimizer(mwe, Gurobi.Optimizer)
set_attribute(mwe, "Method", 2)
set_attribute(mwe, "Crossover", 0)
set_attribute(mwe, "PreDual", 0)


optimize!(mwe)



# User-callback calls 148, time in user-callback 0.00 sec

primal_status(mwe)       # FEASIBLE_POINT::ResultStatusCode = 1
dual_status(mwe)         # FEASIBLE_POINT::ResultStatusCode = 1
termination_status(mwe)  # OPTIMAL::TerminationStatusCode = 1

solution_summary(mwe)

# * Solver : Gurobi

# * Status
#   Result count       : 1
#   Termination status : OPTIMAL
#   Message from the solver:
#   "Model was solved to optimality (subject to tolerances), and an optimal solution is available."

# * Candidate solution (result #1)
#   Primal status      : FEASIBLE_POINT
#   Dual status        : FEASIBLE_POINT
#   Objective value    : 2.76039e+05

# * Work counters
#   Solve time (sec)   : 4.66490e-03
#   Barrier iterations : 30
#   Node count         : 0

has_values(mwe)  # true
has_duals(mwe)   # true

objective_value(mwe)       # 276039.476076866 (true value: 2.7704847864e+05)
dual_objective_value(mwe)  # ERROR: Gurobi Error 10005: Unable to retrieve attribute 'ObjBound'

objective_function(mwe)    # θ
value(first(keys(objective_function(mwe).terms)))           # 277822.2353388086

# way higher than the tolerances?

set_attribute(mwe, "NumericFocus", 3)
set_attribute(mwe, "OptimalityTol", 1e-9)
set_attribute(mwe, "FeasibilityTol", 1e-9)
set_attribute(mwe, "BarConvTol", 1e-16)
optimize!(mwe)

# prints "Sub-optimal termination - objective 2.76103058e+05"
dual_status(mwe)  # UNKNOWN_RESULT_STATUS::ResultStatusCode = 8

# Same result

_m = Model(Gurobi.Optimizer)
@variable(_m, x >= 0)
@objective(_m, Min, x)
optimize!(_m)

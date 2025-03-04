m_main = direct_model(HiGHS.Optimizer())

m_sub = direct_model(HiGHS.Optimizer())
set_attribute(m_sub, "presolve", "off")

# -------------------------------------------

@variable(m_main, my[1:2] >= 0)
@objective(m_main, Min, sum(my))

@variable(m_sub, sx[1:2] >= 0)
@variable(m_sub, sy[1:2])
sc = @constraints(m_sub, begin
    sx[1] <= sy[1]
    sx[2] <= sy[2]
    sum(sx) >= 5.0
end)
@objective(m_sub, Min, 2 * sx[1] + 3 * sx[2])

# -------------------------------------------

optimize!(m_main)
value.(my)

fix.(sy, value.(my); force=true) # value(my)
optimize!(m_sub)

objective_value(m_sub)
dual_objective_value(m_sub)


λ = dual.(FixRef.(sy))
@constraint(m_main, λ' * (fix_value.(sy) .- my) >= dual_objective_value(m_sub))


constrnump = Ref{Int}()
Gurobi.GRBgetconstrbyname(
    m_sub.moi_backend,
    "fixme",
    pointer_from_objref(constrnump)
)

farkasdual = Ref{Float64}()
Gurobi.GRBgetdblattrelement(
    m_sub.moi_backend,
    "FarkasDual",
    2, #constrnump[],
    pointer_from_objref(farkasdual)
)





model = read_from_file("national_scale_$T.mps"; format = JuMP.MOI.FileFormats.FORMAT_MPS)
set_optimizer(model, HiGHS.Optimizer)

set_attribute(model, "presolve", "off")
optimize!(model)

set_attribute(model, "presolve", "on")
optimize!(model)

of = objective_function(model)
@objective(model, Max, -of)


b = unsafe_backend(model)
b.solution
b.options

set_attribute(model, "simplex_strategy", 1)

set_start_values(model, variable_primal_start = (x) -> 0.0)

con = all_constraints(model; include_variable_in_set_constraints = false)[1]
JuMP.set_normalized_rhs(con, JuMP.normalized_rhs(con))


using JuMP
import HiGHS
import Random: seed!

# Random data.
seed!(12345)
N, M = (500, 500)
A = abs.(randn(N, M))
b = randn(N)
c = abs.(randn(M)), abs.(randn(M))
r = randn(M) .* 100

# Model setup.
m = Model(Gurobi.Optimizer)
# set_attribute(m, "simplex_strategy", 4)
@variable(m, r[i] <= x[i = 1:M] <= (r .+ 100)[i])
@constraint(m, A * x .>= b)
@objective(m, Min, c[1]' * x)

# Solve 1.
optimize!(m)                            # := line 42

# Modify.
@objective(m, Min, c[2]' * x)

# Solve 2.
MOI.Utilities.reset_optimizer(m)
optimize!(m)


# with line 42 => 410 iterations
# w/o line 42  => 241 iterations


set_attribute(model_1, "LPWarmStart", 0)

# Solving LP without presolve, or with basis, or unconstrained
# => 43 simplex iterations

# Solving the presolved LP
# => 36 simplex iterations



# none / none => 15 (no presolve)
# no solve / none => 10 (presolve)

model = direct_model(HiGHS.Optimizer())
@variable(model, x >= 0)

@constraint(model, c1, x <= 1.0)
@constraint(model, c2, x >= 3.0)

set_attribute(model, "presolve", "off")
optimize!(model)

dual_status(model)  # INFEASIBILITY_CERTIFICATE::ResultStatusCode = 4

set_normalized_rhs(c1, 2.9999)
set_attribute(model, "presolve", "on")
optimize!(model)

primal_status(model)
dual_status(model)
termination_status(model)














using JuMP
import Gurobi

mwe = read_from_file("mwe.mps"; format = MOI.FileFormats.FORMAT_MPS)
set_optimizer(mwe, Gurobi.Optimizer)
set_attribute(mwe, "Method", 2)
set_attribute(mwe, "Crossover", 0)

set_attribute(mwe, "NumericFocus", 3)

set_attribute(mwe, "PreDual", 0)
set_attribute(mwe, "BarConvTol", 1e-2)

optimize!(mwe)

MOI.get(mwe, Gurobi.ModelAttribute("ObjBound"))


objective_value(mwe) == value(objective_function(mwe))

# From iteration 29
ub = 2.77822235e+05     # value(objective_function(mwe))
lb = 2.76039476e+05     # objective_value(mwe)
gap_abs = ub - lb       # ~ 1782.759
gap_rel = gap_abs / ub  # ~ 0.006417

MOI.get(mwe, MOI.RawOptimizerAttribute("BarConvTol"))

objective_value(mwe)

MOI.get(mwe, Gurobi.ModelAttribute("Status"))  # = 2 ("OPTIMAL")

# Why does that not return 13 ("SUBOPTIMAL"), [which states](https://docs.gurobi.com/projects/optimizer/en/current/reference/numericcodes/statuscodes.html#secstatuscodes) _"Unable to satisfy optimality tolerances; a sub-optimal solution is available."_
[BarConvTol](https://www.gurobi.com/documentation/current/refman/barconvtol.html) states _"The barrier solver terminates when the relative difference between the primal and dual objective values is less than the specified tolerance (with a GRB_OPTIMAL status)."_ The gap I see in the log is greater than 6.4e-3, while the tolerance defaults to 1e-8.

The dual violations, as you showed, are still quite large - shouldn't [OptimalityTol](https://www.gurobi.com/documentation/current/refman/optimalitytol.html) also prevent flagging it as optimal?

MOI.get(mwe, Gurobi.ModelAttribute("ObjVal"))  # 276039.476076866
MOI.get(mwe, Gurobi.ModelAttribute("ObjBound"))  # Gurobi Error 10005



Can you confirm that my guess is correct, that with default settings (`PreDual = -1`) it chooses to solve the dual?

I can only expect the true objective to be between the lower and upper bound. With the problem being a minimization problem, if I'm looking at a statement of "we know an optimal decision will at most cost us X", then I need access to the upper bound, which I'd get from the primal objective value (if it is actually solving the dual in the background, this then would be the dual objective value of what is being solved). Even though some intermediate value is available (as evident in the log), it's not accessible, because I can't access `ObjBound`.



Test case `PreDual=0`

```julia
set_attribute(mwe, "PreDual", 0)
set_attribute(mwe, "BarConvTol", 1e-2)

# Last iteration in log:
#   13   2.76998764e+05  2.76994994e+05  0.00e+00 2.97e-08  2.19e-02     0s

MOI.get(mwe, Gurobi.ModelAttribute("ObjVal"))  # 276998.7644756779
objective_value(mwe)                           # 276998.7644756779
value(objective_function(mwe))                 # 276998.7644756779
```

Test case `PreDual=1`

```julia
set_attribute(mwe, "PreDual", 1)
set_attribute(mwe, "BarConvTol", 1e-2)

# Last iteration in log:
#   14   2.75641249e+05  2.77856133e+05  4.87e-06 5.13e-06  3.22e+00     0s

MOI.get(mwe, Gurobi.ModelAttribute("ObjVal"))  # 275641.248734008
objective_value(mwe)                           # 275641.248734008
value(objective_function(mwe))                 # 277856.1332178363
```







using JuMP
import HiGHS
import Random: seed!

# Random data.
seed!(12345)
N, M = (500, 500)
A = abs.(randn(N, M))
b = randn(N)
c = abs.(randn(M)), abs.(randn(M))
r = randn(M) .* 100

# Model setup.
m = Model(Gurobi.Optimizer)
# set_attribute(m, "BarConvTol", 1e-4)

# set_attribute(m, "simplex_strategy", 4)
@variable(m, r[i] <= x[i = 1:M] <= (r .+ 100)[i])
@constraint(m, A * x .>= b)
@objective(m, Min, c[1]' * x + x[1] * x[1])

# Solve 1.
optimize!(m)

objective_bound(m)

objective_value(m)
dual_objective_value(m)









# TODO TOMORROW: post & only use primal for HiGHS for now ...


using JuMP
import HiGHS

mwe = read_from_file("mwe_highs_optimal_and_infeasible.mps")
set_optimizer(mwe, HiGHS.Optimizer)
set_attribute(mwe, "solver", "ipm")
set_attribute(mwe, "presolve", "off")
set_attribute(mwe, "run_crossover", "off")
set_attribute(mwe, "dual_feasibility_tolerance", 1e-10)

optimize!(mwe)

p = Ref{HiGHS.HighsInt}()
HiGHS.Highs_getIntInfoValue(unsafe_backend(mwe), "crossover_iteration_count", p)
p

HiGHS.Highs_getIntInfoValue

senseP = Ref{Int32}()
ret = HiGHS.Highs_getObjectiveSense(unsafe_backend(mwe), senseP)

set_attribute(mwe, "run_crossover", "off")
set_attribute(mwe, "dual_feasibility_tolerance", 1e-10)

mwe = subs(model)[end-1]


dual.(all_constraints(mwe; include_variable_in_set_constraints = true))


termination_status(mwe)

primal_status(mwe)
dual_status(mwe)

write_to_file(mwe, "mwe_highs_optimal_and_infeasible.mps")

# HiGHS log output
# Running HiGHS 1.7.2 (git hash: 5ce7a2753): Copyright (c) 2024 HiGHS under MIT licence terms
# Coefficient ranges:
#   Matrix [1e-09, 1e+07]
#   Cost   [8e-02, 2e+00]
#   Bound  [1e+03, 1e+06]
#   RHS    [1e-13, 4e+11]
# WARNING: Problem has excessively large bounds: consider scaling the bounds by 1e+2 or less, or setting option user_bound_scale to -6 or less
# Presolving model
# 4071 rows, 179 cols, 23438 nonzeros  0s
# 875 rows, 179 cols, 8004 nonzeros  0s
# 853 rows, 157 cols, 7960 nonzeros  0s
# Presolve : Reductions: rows 853(-3350); columns 157(-37); elements 7960(-15816)
# Solving the presolved LP
# IPX model has 853 rows, 157 columns and 7960 nonzeros
# Input
#     Number of variables:                                157
#     Number of free variables:                           0
#     Number of constraints:                              853
#     Number of equality constraints:                     0
#     Number of matrix entries:                           7960
#     Matrix range:                                       [1e-09, 1e+07]
#     RHS range:                                          [1e-13, 4e+11]
#     Objective range:                                    [4e-02, 2e+00]
#     Bounds range:                                       [1e+03, 1e+06]
# Preprocessing
#     Dualized model:                                     yes
#     Number of dense columns:                            0
#     Range of scaling factors:                           [4.88e-04, 5.12e+02]
# IPX version 1.0
# Interior Point Solve
#  Iter     P.res    D.res            P.obj           D.obj        mu     Time
#    0   1.61e+00 4.06e+08  -3.62747966e+06 -2.73315498e+05  7.27e+08       0s
#  Constructing starting basis...
#    1   7.45e-01 1.49e+08   2.16503232e+10 -5.82216633e+09  3.36e+08       0s
#    2   3.35e-01 1.07e+08   2.63366163e+10 -6.81053553e+09  2.06e+08       0s
#    3   1.51e-01 6.50e+07   2.48884472e+10 -7.62696235e+09  1.26e+08       0s
#    4   5.94e-02 4.45e+07   1.81062320e+10 -7.24209359e+09  7.43e+07       0s
#    5   1.07e-02 6.28e+06   9.16378252e+09 -3.88644637e+09  1.79e+07       0s
#    6   9.24e-04 7.96e+05   1.75717528e+09 -7.63840536e+08  2.64e+06       0s
#    7   1.46e-05 6.36e+04   1.34956075e+08 -7.48031367e+07  1.92e+05       0s
#    8   6.84e-06 9.10e+03   6.44043972e+07 -1.30155029e+07  6.90e+04       0s
#    9   9.01e-07 9.10e-03   8.93610184e+06 -2.25615838e+06  9.87e+03       0s
#   10   4.00e-07 6.86e-03   6.09280309e+06 -2.24150381e+06  7.32e+03       0s
#   11   1.23e-07 3.12e-04   1.90785320e+06 -1.07336378e+06  2.61e+03       0s
#   12   3.72e-08 7.05e-05   3.65024020e+05 -6.48372243e+05  8.88e+02       0s
#   13   4.57e-09 4.04e-05   1.77437232e+04 -5.72094405e+05  5.16e+02       0s
#   14   5.51e-10 2.40e-05  -1.64894005e+05 -4.65095112e+05  2.62e+02       0s
#   15   1.78e-10 8.46e-06  -2.26382751e+05 -3.39900358e+05  9.92e+01       0s
#   16   1.37e-10 1.67e-06  -2.34368908e+05 -3.07124631e+05  6.36e+01       0s
#   17   4.64e-11 4.17e-07  -2.58913594e+05 -2.89300544e+05  2.66e+01       0s
#   18   4.96e-12 2.38e-07  -2.68725857e+05 -2.83450094e+05  1.29e+01       0s
#   19   9.42e-13 2.38e-07  -2.74756709e+05 -2.79008239e+05  3.72e+00       0s
#   20   6.36e-13 1.19e-07  -2.75326767e+05 -2.77990200e+05  2.33e+00       0s
#   21   4.23e-13 1.19e-07  -2.75827237e+05 -2.77402153e+05  1.38e+00       0s
#   22   8.39e-14 1.79e-07  -2.76758400e+05 -2.77139768e+05  3.33e-01       0s
#   23   3.73e-14 1.79e-07  -2.76902232e+05 -2.77081895e+05  1.57e-01       0s
#   24   7.77e-15 1.79e-07  -2.77007716e+05 -2.77061901e+05  4.74e-02       0s
#   25   1.22e-15 2.98e-07  -2.77037652e+05 -2.77051172e+05  1.18e-02       0s
#   26   5.64e-16 2.38e-07  -2.77045541e+05 -2.77048912e+05  2.95e-03       0s
#   27   3.89e-16 1.79e-07  -2.77047733e+05 -2.77048357e+05  5.45e-04       0s
#   28   3.21e-16 2.38e-07  -2.77048073e+05 -2.77048209e+05  1.19e-04       0s
#   29   4.09e-16 1.19e-07  -2.77048176e+05 -2.77048192e+05  1.39e-05       0s
#   30*  4.41e-16 1.79e-07  -2.77048188e+05 -2.77048188e+05  3.20e-07       0s
# Summary
#     Runtime:                                            0.02s
#     Status interior point solve:                        optimal
#     Status crossover:                                   not run
#     objective value:                                    2.77048188e+05
#     interior solution primal residual (abs/rel):        1.83e-04 / 5.17e-16
#     interior solution dual residual (abs/rel):          8.86e-10 / 3.53e-10
#     interior solution objective gap (abs/rel):          3.76e-04 / 1.36e-09
# Ipx: IPM       optimal
# Model   status      : Optimal
# IPM       iterations: 30
# Objective value     :  2.7704818697e+05
# HiGHS run time      :          0.02

# Solution summary
# * Solver : HiGHS

# * Status
#   Result count       : 1
#   Termination status : OPTIMAL
#   Message from the solver:
#   "kHighsModelStatusOptimal"

# * Candidate solution (result #1)
#   Primal status      : FEASIBLE_POINT
#   Dual status        : FEASIBLE_POINT
#   Objective value    : 2.77048e+05
#   Objective bound    : 0.00000e+00
#   Relative gap       : Inf
#   Dual objective value : -2.22082e+05

# * Work counters
#   Solve time (sec)   : 2.48494e-02
#   Simplex iterations : 0
#   Barrier iterations : 30
#   Node count         : -1

unsafe_backend(mwe).solution.status

unsafe_backend(mwe).solution.model_status

# nope
unsafe_backend(mwe).solution.dual_solution_status   # kSolutionStatusFeasible
unsafe_backend(mwe).solution.primal_solution_status   # kSolutionStatusFeasible


MOI.get(mwe, MOI.DualObjectiveValue(1))

HiGHS.Highs_getObjectiveValue(unsafe_backend(mwe))

objective_bound(mwe)

p = Ref{Cdouble}()
HiGHS.Highs_getDoubleInfoValue(unsafe_backend(m_main), "mip_dual_bound", p)
p


using JuMP
import HiGHS
import Random: seed!

# Random data.
seed!(12345)
N, M = (500, 500)
A = abs.(randn(N, M))
b = randn(N)
c = abs.(randn(M)), abs.(randn(M))
r = randn(M) .* 100

# Model setup.
m = Model(HiGHS.Optimizer)
set_attribute(m, "solver", "ipm")
set_attribute(m, "run_crossover", "off")

@variable(m, r[i] <= x[i = 1:M] <= (r .+ 100)[i])
@constraint(m, A * x .>= b)
@objective(m, Min, c[1]' * x + x[1] * x[1])

optimize!(m)

MOI.get(m, MOI.DualObjectiveValue(1))





import HiGHS
model = Model(HiGHS.Optimizer);
@variable(model, x >= 1);
@objective(model, Min, 2 * x + 1);
set_attribute(mwe, "solver", "ipm")
set_attribute(mwe, "run_crossover", "off")
optimize!(model)
dual_objective_value(m)

p = Ref{Cdouble}()
HiGHS.Highs_getDoubleInfoValue(unsafe_backend(model), "mip_dual_bound", p)
p




con = all_constraints(mwe; include_variable_in_set_constraints = false)

π = dual.(con)
b = normalized_rhs.(con)

π' * b


F = MOI.get(mwe, MOI.ObjectiveFunctionType())
f = MOI.get(mwe, MOI.ObjectiveFunction{F}())
MOI.constant(f, Float64)
MOI.get(mwe, MOI.ObjectiveSense())
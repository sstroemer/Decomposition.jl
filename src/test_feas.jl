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
m = Model(HiGHS.Optimizer)
set_attribute(m, "simplex_strategy", 4)
@variable(m, r[i] <= x[i = 1:M] <= (r .+ 100)[i])
@constraint(m, A * x .>= b)
@objective(m, Min, c[1]' * x)

# Solve 1.
optimize!(m)                            # := line 42

# Modify.
@objective(m, Min, c[2]' * x)

# Solve 2.
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

using JuMP
import HiGHS, Gurobi
using Random: seed!

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seed!(1234)

N = 100
A = rand(N, N) .- 0.5
b = rand(N)
c = rand(N)

M = 10
_A = zeros(M, N)
for i in 1:M
    _A[i, i] = -1
end
A = vcat(A, _A)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function sub_primal(model, x)
    @variable(model, 0 <= y[1:N])
    link = @constraint(model, A * y .>= vcat(b, -x))
    @objective(model, Min, c' * y)
    optimize!(model)

    return dual_objective_value(model), dual.(link[(N+1):end])
end

function sub_dual(model, x)
    @variable(model, 0 <= y[1:(N+M)])
    @constraint(model, A' * y .<= c)
    @objective(model, Max, vcat(b, -x)' * y)
    optimize!(model)

    return objective_value(model), value.(y[(N+1):end])
end

# Gurobi: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#infunbdinfo
# When this parameter is set additional information will be computed when a model is determined to be infeasible or
# unbounded, and a simplex basis is available (from simplex or crossover). Note that if a model is determined to be
# infeasible or unbounded when solving with barrier, prior to crossover, then this additional information will not be
# available.

# => A basic solution is necessary for both primal and dual rays.
# => Barrier cannot be used, even when crossover is enabled.

opt = Gurobi.Optimizer

main = Model(opt)
@variable(main, x[1:M] >= 0)
@objective(main, Min, sum(x))

k = 0
for i in 1:100
    optimize!(main)
    vx = value.(x)

    sub = Model(opt)

    # Note on HiGHS IPM and crossover: Crossover will often times be skipped. However, it will sometimes encounter
    # `WARNING: IPX solution is imprecise, so clean up with simplex`. Without crossover it would in these cases fail.

    set_attributes(sub, "InfUnbdInfo" => 1, "DualReductions" => 0, "Method" => 0)            # primal simplex
    # set_attributes(sub, "InfUnbdInfo" => 1, "DualReductions" => 0, "Method" => 1)            # dual simplex
    # set_attributes(sub, "presolve" => "off", "solver" => "simplex", "simplex_strategy" => 4) # primal simplex
    # set_attributes(sub, "presolve" => "off", "solver" => "simplex", "simplex_strategy" => 1) # dual simplex
    # set_attributes(sub, "presolve" => "off", "solver" => "ipm", "run_crossover" => "on")

    obj, λ = sub_primal(sub, vx)
    # obj, λ = sub_dual(sub, vx)

    if is_solved_and_feasible(sub)
        break
    end

    @constraint(main, obj - λ' * (x - vx) <= 0)
    k += 1
end

# Number of iterations:
# Gurobi, solving primal, primal simplex: 26
# Gurobi, solving primal, dual simplex  : 21
# Gurobi, solving dual,   primal simplex: 26
# Gurobi, solving dual,   dual simplex  : 30
# HiGHS,  solving primal, primal simplex: 39
# HiGHS,  solving primal, dual simplex  : 20
# HiGHS,  solving primal, ipm           : 44
# HiGHS,  solving dual,   primal simplex: 42
# HiGHS,  solving dual,   dual simplex  : 44
# HiGHS,  solving dual,   ipm           : 44

_presolve_papilo(node::AbstractNode) = _presolve_papilo(to_file(node))

function _presolve_papilo(node::AbstractFileNode)
    @error "Presolve is currently deactivated, since it is 'wonky' when interested in duals (and not resulting in a good presolve then), and otherwise only working for primal results: see https://github.com/scipopt/papilo/issues/50"
    return node
    # SCIP_PaPILO_jll.papilo() do exe
    #     run(`$exe solve -f tmp/test.mps -v tmp/test.postsolve -r tmp/test_presolved.mps -x tmp/soplex.set --parameter-settings tmp/papilo.set `)
    # end

    # node = "tmp/test_presolved.mps" |> to_model
    # set_optimizer(node.model, HiGHS.Optimizer)
    # optimize!(node.model)
    # Decomposition._write_papilo_sol(node.model, "tmp/test_presolved")

    # SCIP_PaPILO_jll.papilo() do exe
    #     run(`$exe postsolve -v tmp/test.postsolve -u tmp/test_presolved.sol --dual-reduced-solution tmp/test_presolved.dual.sol --costs-reduced-solution tmp/test_presolved.reducedcost.sol -l tmp/test.sol --dualsolution tmp/test.dual.sol -c tmp/test.reducedcost.sol`)
    # end

    @debug "_presolve_papilo(::AbstractNode)" node

    filename_full = node.filename
    raw_filename = rsplit(filename_full, "."; limit=2)[1]
    filename_reduced = "$(raw_filename).papilo.reduced.mps"
    filename_postsolve = "$(raw_filename).papilo.postsolve"

    # These are prepared for the postsolve step.
    filename_reduced_sol = "$(raw_filename).papilo.reduced.sol"
    filename_full_sol = "$(raw_filename).papilo.full.sol"

    Suppressor.@suppress PaPILO.presolve_write_from_file(filename_full, filename_postsolve, filename_reduced)

    return Node{FileNodePresolved}(node; filename=filename_reduced, filename_original=filename_full, filename_postsolve=filename_postsolve, filename_reduced_sol=filename_reduced_sol, filename_full_sol=filename_full_sol)
end

function _postsolve_papilo(node::FileNodePresolved)
    SCIP_PaPILO_jll.papilo() do exe
        run(`$(exe) postsolve -v $(node.filename_postsolve) -u $(node.filename_reduced_sol) -l $(node.filename_full_sol)`)
    end
end

function _write_papilo_sol(model::JuMP.Model, filename::String)
    all_vars = JuMP.all_variables(model)
    all_constr = JuMP.all_constraints(model; include_variable_in_set_constraints=false)
    open("$(filename).sol", "w") do solfile
        for var in all_vars
            write(solfile, "$(JuMP.name(var))\t$(JuMP.value(var))\n")
        end
    end
    open("$(filename).reducedcost.sol", "w") do solfile
        for var in all_vars
            write(solfile, "$(JuMP.name(var))\t$(JuMP.reduced_cost(var))\n")
        end
    end
    open("$(filename).dual.sol", "w") do solfile
        for con in all_constr
            write(solfile, "$(JuMP.name(con))\t$(JuMP.dual(con))\n")
        end
    end
end

function _resolve_ranged_constraints!(model::JuMP.Model)
    ranged_constraints = JuMP.all_constraints(model, JuMP.VariableRef, MOI.Interval{Float64})
    isempty(ranged_constraints) || @warn "Converting ranged constraints to bounds, to 'fix' dualization"

    for rc in ranged_constraints
        func = JuMP.constraint_object(rc).func
        !(func isa JuMP.VariableRef) && @critical "Ranged constraint resolution currenlty only implemented for variable bounds"
        set = JuMP.constraint_object(rc).set
        JuMP.delete(model, rc)  # TODO: precalc the upper/lower and do the removal in one pass beforehand
        JuMP.set_lower_bound(func, set.lower)
        JuMP.set_upper_bound(func, set.upper)
    end    

    return nothing
end

function _to_program_lp(node::AbstractModelNode)
    @debug "_to_program_lp(::AbstractModelNode)" node

    @warn "_to_program_lp does not properly handle presolved models (probably ranged constraints)"

    model = node.model
    lpmd = JuMP.lp_matrix_data(model)

    # Get objective function.
    c = (lpmd.sense == MOI.MIN_SENSE) ? lpmd.c : -lpmd.c
    d = lpmd.c_offset

    # Get sign of constraints to split into `Ax ≤ b` and `Gx = h`.
    idx_eq = (lpmd.b_lower .== lpmd.b_upper)
    idx_leq = (.!idx_eq) .& (lpmd.b_upper .< Inf)
    idx_geq = (.!idx_eq) .& (lpmd.b_lower .> -Inf)

    A = vcat(lpmd.A[idx_leq, :], -lpmd.A[idx_geq, :])
    b = vcat(lpmd.b_upper[idx_leq], -lpmd.b_lower[idx_geq])

    G = lpmd.A[idx_eq, :]
    h = lpmd.b_upper[idx_eq]

    # Prepare non-negativity.
    lb_x = lpmd.x_lower

    # Get all variables `x <= ub`, and replace them by `y = -x`, with `y >= -ub`.
    # That means reversing the sign in `c`, `A`, and `G`.
    ub_x = lpmd.x_upper
    reverse_vars = ub_x .!= Inf
    c[reverse_vars] = -c[reverse_vars]
    A[:, reverse_vars] = -A[:, reverse_vars]
    G[:, reverse_vars] = -G[:, reverse_vars]
    ub_x = -ub_x

    # Merge the true lower bounds and the reversed upper bounds.
    lb = max.(lb_x, ub_x)

    # Ensure that `lb` only concerns non-negativity.
    idx_bound_to_constraint = (lb .!= 0) .& (lb .!= -Inf)
    if any(idx_bound_to_constraint)
        @error "Non-zero variable bound to constraint conversion not implemented"
    end

    # TODO: retain names
    return Node{ProgramNodeLP}(
        node;
        program = ProgramLP(A, G, b, c, d, h, lb, size(A, 1), size(G, 1), size(A, 1) + size(G, 1), size(A, 2)),
    )
end

to_program(node::String) = to_program(from_file(node))
to_program(node::AbstractNode) = to_program(to_model(node))
to_program(node::AbstractProgramNode) = node
function to_program(node::AbstractModelNode)
    @debug "to_program(::AbstractModelNode)" node

    return _to_program_lp(node)

    # TODO: _is_lp does not account for: range constraints (ok), QPs (not ok), QCPs (not ok), ...
    #       Code as it was planned:
    # 
    # JuMP._is_lp(node.model) && (return _to_program_lp(node))
    #
    # @critical "`to_program` not implemented for model type, currently supported: LP"
end

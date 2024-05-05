"""
Linear Program (LP) with `n` variables, and `m` constraints, in the form:

    min c'x + d
    s.t.
        Ax ≤ b
        Gx = h
        x ≥ lb,     lb_i ∈ {-∞, 0}, ∀ i ∈ M

Where:
    x ∈ R^n
    
    c ∈ R^n
    d ∈ R
    
    A ∈ R^{m1 x n}
    b ∈ R^m1

    G ∈ R^{m2 x n}
    h ∈ R^m2
    
    M := {1, ..., m}
    N := {1, ..., n}

With:
    m := m1 + m2
"""
struct ProgramLP <: AbstractProgram
    A::SparseArrays.SparseMatrixCSC{Float64, Int64}
    G::SparseArrays.SparseMatrixCSC{Float64, Int64}

    b::SparseArrays.SparseVector{Float64, Int64}
    c::SparseArrays.SparseVector{Float64, Int64}
    d::Float64
    h::SparseArrays.SparseVector{Float64, Int64}

    lb::SparseArrays.SparseVector{Float64, Int64}

    m1::Int
    m2::Int
    m::Int
    n::Int
end

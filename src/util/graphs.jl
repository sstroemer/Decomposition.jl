"""
Helper function to get the graph structure.
"""
function create_adjacency_matrix(A::SparseArrays.SparseMatrixCSC, rm::Set{Int64})
    tA = SparseArrays.SparseMatrixCSC((A .!= 0)')
    rvs = SparseArrays.rowvals(tA)
    
    I = Int64[]
    J = Int64[]

    for i in 1:size(A, 1)
        nodes = rvs[SparseArrays.nzrange(tA, i)]
        for j in eachindex(nodes)
            u = nodes[j]
            (u in rm) && continue
            for k in (j+1):length(nodes)
                v = nodes[k]
                (v in rm) && continue
                push!(I, u)
                push!(J, v)
            end
        end
    end

    return LinearAlgebra.Symmetric(SparseArrays.sparse(I, J, 1, size(A, 2), size(A, 2)))
end

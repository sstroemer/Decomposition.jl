_presolve_papilo(node::AbstractNode) = _presolve_papilo(to_file(node))

function _presolve_papilo(node::AbstractFileNode)
    @debug "_presolve_papilo(::AbstractNode; ::Bool)" node

    filename_full = node.filename
    raw_filename = rsplit(filename_full, "."; limit=2)[1]
    filename_reduced = "$(raw_filename).papilo.reduced.mps"
    filename_postsolve = "$(raw_filename).papilo.postsolve"

    PaPILO.presolve_write_from_file(filename_full, filename_postsolve, filename_reduced)

    return Node{NodeModelPresolved}(node; model=JuMP.read_from_file(filename_reduced), filename_full=filename_full, filename_reduced=filename_reduced, filename_postsolve=filename_postsolve)
end

# Decomposition.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sstroemer.github.io/Decomposition.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sstroemer.github.io/Decomposition.jl/dev/)
[![Build Status](https://github.com/sstroemer/Decomposition.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sstroemer/Decomposition.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sstroemer/Decomposition.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sstroemer/Decomposition.jl)

General purpose package to deal with various decomposition approaches for optimization problems.

## Architecture
TBA

## Visualizing a workflow

A workflow, e.g. created using `Workflow(node)` on any node, can be printed into a tree structure using `print_workflow(workflow)`:

```julia
julia> print_workflow(Workflow(node))

Root [00000]
╤═══════════
FileGeneral [00001]: input_file.mps
╰─ ModelGeneral [00002]: 5206 vars and 7854 cons
   ╰─ FileGeneral [00003]: tmpfile.mps
      ╰─ ModelGeneral [00004]: 5206 vars and 7854 cons
         ╰─ SimpleBendersDecomposition [00005]
            ├─ ModelGeneral [00006]: 5206 vars and 7854 cons
            ├─ ModelGeneral [00007]: 5206 vars and 7854 cons
            ┆  ╰─ SimpleBendersDecomposition [0000b]
            ┆     ├─ ModelGeneral [0000c]: 5206 vars and 7854 cons
            ┆     ├─ ModelGeneral [0000d]: 5206 vars and 7854 cons
            ┆     ├─ ModelGeneral [0000e]: 5206 vars and 7854 cons
            ┆     ├─ ModelGeneral [0000f]: 5206 vars and 7854 cons
            ┆     ╰─ ModelGeneral [00010]: 5206 vars and 7854 cons
            ├─ ModelGeneral [00008]: 5206 vars and 7854 cons
            ├─ ModelGeneral [00009]: 5206 vars and 7854 cons
            ╰─ ModelGeneral [0000a]: 5206 vars and 7854 cons
```

## Details on Nodes

The following high-level node types currently exist:
- `AbstractNode`: Abstract base type for all nodes
- `AbstractFileNode`: Abstract base type for nodes that read/write files
- `AbstractModelNode`: Abstract base type for nodes that work with `JuMP` models
- `AbstractProgramNode`: Abstract base type for nodes that work with mathematical programs (e.g., LPs)

More details on the specific types are given below.

### `NodeRoot`
TBA

### `FileNodeGeneral`
TBA

### `FileNodePresolved`
TBA

### `ModelNodeGeneral`
TBA

### `ModelNodeDualization`

A `ModelNodeDualization` node relies on an initial `JuMP` model from its parent node; if that is not satisfied, `Decomposition.jl` will try to create one. Based on that given model, it then uses a `MockOptimizer` (from `MOI.Utilities`), and wraps it into a `DualOptimizer` using `Dualization.jl`. On triggering a solve, it will either:

- replace the `MockOptimizer` with a proper one, if the node is a final leaf node, or
- use the `MockOptimizer` linked to results from the node's child

This allows extracting primal/dual results from the "original model" (= in the way it was specified by the node's parent) based on the built-in conversion functionality in `Dualization.jl`.

### `ProgramNodeLP`
TBA

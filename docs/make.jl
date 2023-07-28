using Decomposition
using Documenter

DocMeta.setdocmeta!(Decomposition, :DocTestSetup, :(using Decomposition); recursive=true)

makedocs(;
    modules=[Decomposition],
    authors="sstroemer <stefan.stroemer@gmail.com>",
    repo="https://github.com/sstroemer/Decomposition.jl/blob/{commit}{path}#{line}",
    sitename="Decomposition.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sstroemer.github.io/Decomposition.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sstroemer/Decomposition.jl",
    devbranch="main",
)

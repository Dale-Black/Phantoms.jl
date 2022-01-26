using Phantoms
using Documenter

DocMeta.setdocmeta!(Phantoms, :DocTestSetup, :(using Phantoms); recursive=true)

makedocs(;
    modules=[Phantoms],
    authors="Dale-Black <djblack@uci.edu> and contributors",
    repo="https://github.com/Dale-Black/Phantoms.jl/blob/{commit}{path}#{line}",
    sitename="Phantoms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dale-Black.github.io/Phantoms.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/Dale-Black/Phantoms.jl", devbranch="main")

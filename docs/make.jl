using SymTridiagonalQL
using Documenter

DocMeta.setdocmeta!(SymTridiagonalQL, :DocTestSetup, :(using SymTridiagonalQL); recursive=true)

makedocs(;
    modules=[SymTridiagonalQL],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/jagot/SymTridiagonalQL.jl/blob/{commit}{path}#{line}",
    sitename="SymTridiagonalQL.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/SymTridiagonalQL.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jagot/SymTridiagonalQL.jl",
    devbranch="main",
)

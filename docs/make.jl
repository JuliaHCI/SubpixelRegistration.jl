using SubpixelRegistration
using Documenter

DocMeta.setdocmeta!(SubpixelRegistration, :DocTestSetup, :(using SubpixelRegistration); recursive=true)

makedocs(;
    modules=[SubpixelRegistration],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaHCI/SubpixelRegistration.jl/blob/{commit}{path}#{line}",
    sitename="SubpixelRegistration.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaHCI.github.io/SubpixelRegistration.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaHCI/SubpixelRegistration.jl",
    devbranch="main",
)

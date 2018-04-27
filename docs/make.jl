using Documenter, SubpixelRegistration

makedocs(modules= subpixelRegistration)
deploydocs(
    julia = "nightly",
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/romainFr/SubpixelRegistration.jl.git"
)

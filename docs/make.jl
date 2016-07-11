using Documenter, subpixelRegistration

makedocs(modules= subpixelRegistration)
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/romainFr/subpixelRegistration.jl.git"
)

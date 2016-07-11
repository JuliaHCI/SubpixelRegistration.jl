using Documenter, subpixelRegistration

makedocs()
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/romainFr/subpixelRegistration.jl.git"
)

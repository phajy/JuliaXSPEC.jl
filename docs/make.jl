using Documenter
using JuliaXSPEC

makedocs(
    modules = [JuliaXSPEC],
    sitename = "JuliaXSPEC.jl",
    authors = "Andy Young",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://phajy.github.io/JuliaXSPEC.jl",
        mathengine = Documenter.MathJax3(),
    ),
    pages = [
        "Home" => "index.md",
        "Installation and build" => "build.md",
        "Quick start" => "quickstart.md",
        "Concepts" => [
            "Units, bins and normalisation" => "units.md",
            "Convolution" => "convolution.md",
            "Caching on a grid" => "caching.md",
        ],
        "Writing a model" => "models.md",
        "Verification" => "verification.md",
        "Reference" => "reference.md",
    ],
    checkdocs = :exports,
    warnonly = [:cross_references],
)

deploydocs(
    repo = "github.com/phajy/JuliaXSPEC.jl.git",
    push_preview = true,
)

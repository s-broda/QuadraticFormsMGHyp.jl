using Documenter, QuadraticFormsMGHyp, DocThemeIndigo
indigo = DocThemeIndigo.install(QuadraticFormsMGHyp)
makedocs(;
    modules=[QuadraticFormsMGHyp],
    format=Documenter.HTML(assets=String[indigo]),
    pages=[
        "Home" => "index.md",
    ],
    sitename="QuadraticFormsMGHyp.jl",
    authors="Simon Broda",
)

deploydocs(;
    repo="github.com/s-broda/QuadraticFormsMGHyp.jl",
)

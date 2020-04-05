using Documenter, QuadraticFormsMGHyp

makedocs(;
    modules=[QuadraticFormsMGHyp],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/s-broda/QuadraticFormsMGHyp.jl/blob/{commit}{path}#L{line}",
    sitename="QuadraticFormsMGHyp.jl",
    authors="Simon Broda",
)

deploydocs(;
    repo="github.com/s-broda/QuadraticFormsMGHyp.jl",
)

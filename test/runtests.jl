using QuadraticFormsMGHyp
import QuadraticFormsMGHyp: lklam
using Test
using LinearAlgebra
using FinancialToolbox
using ToeplitzMatrices
using Random
using Plots

@testset "2sls" begin
    include("2sls.jl")
    spacdf[isnan.(spacdf)] .= 0
    @test sum(abs2.(cdf - spacdf)) / length(cdf) < 0.001
end
@testset "es" begin
    include("es.jl")
    @test sum(abs2.(ccdf - ccdfspa)) / length(ccdf) < 0.0001
    @test sum(abs2.(pm - pmspa)) / length(pm) < 0.1
end

@static if ~haskey(ENV, "CI")
    pgfplotsx()
    spacdf[abs.(cdf .- spacdf) .> .1] .= NaN # filter errors due to nonsingularity around the mean
    p1 = plot(xvec, cdf', color=collect(1:length(nuvec))', lw=1.25, legend=:topleft, labels=permutedims("\$\\nu =".*string.(nuvec).*"\$"))
    plot!(xvec, spacdf', color=collect(1:length(nuvec))', lw=1.25, labels="", ls=:dash)
    xlabel!("\$\\hat{\\beta}_{2SLS}\$")
    ylabel!("\$f(\\hat{\\beta}_{2SLS})\$")
    title!("Exact (solid) and approximate (dashes) densities")

    p2 = plot(ccdf, pm, labels="", lw=1.25)
    plot!(ccdfspa, pmspa, ls=:dash, lw=1.25, labels="")
    xlabel!("VaR level \$\\alpha\$")
    ylabel!("\$ES_L^{(\\alpha)}\$")
    title!("Exact (solid) and approximate (dashes) expected shortfall")

    plot(p1, p2, layout=@layout[a{1.25w, 1.25h} b{1.25w, 1.25h}])
    savefig("figure1.pdf")
end

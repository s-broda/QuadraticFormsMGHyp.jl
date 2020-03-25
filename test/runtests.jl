using QuadraticFormsMGHyp
using Test
using Plots
using LinearAlgebra
using FinancialToolbox
using ToeplitzMatrices
using Random
pgfplotsx()
if haskey(ENV, "CI")
    doplot = false
else
    doplot = true
end
@testset "2sls" begin
    include("2sls.jl")
    @test isapprox(spapdf[1, 1], 0.017712289301918815, atol=1e-4)
    @test isapprox(pdf[1, 1], 0.01713868502672966, rtol=1e-4)
end
@testset "es" begin
    include("es.jl")
    @test isapprox(ccdf[1, 1], 0.10068491359055437, rtol=1e-4)
    @test isapprox(pm[1, 1], 5.724342391411737, rtol=1e-4)
    @test isapprox(ccdfspa[1, 1], 0.10068589716980002, rtol=1e-4)
    @test isapprox(pmspa[1,1 ], 5.56146456788597, rtol=1e-4)
end

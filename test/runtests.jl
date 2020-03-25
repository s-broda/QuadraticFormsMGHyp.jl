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
    @test spapdf[1, 1] ≈ 0.017712289301918815
    @test pdf[1, 1] ≈ 0.01713868502672966
end
@testset "es" begin
    include("es.jl")
    @test ccdf[1, 1] ≈ 0.10068491359055437
    @test pm[1, 1] ≈ 5.724342391411737
    @test ccdfspa[1, 1] ≈ 0.10068589716980002
    @test pmspa[1,1 ] ≈ 5.56146456788597
end

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
    @test sum(abs2.(cdf - spacdf)) / length(cdf) < 0.001    
end
@testset "es" begin
    include("es.jl")
    @test sum(abs2.(ccdf - ccdfspa)) / length(ccdf) < 0.0001
    @test sum(abs2.(pm - pmspa)) / length(pm) < 0.1
end

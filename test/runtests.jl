using QuadraticFormsMGHyp
using Test
using Plots
using LinearAlgebra
using ToeplitzMatrices

# #@testset "QuadraticFormsMGHyp.jl" begin
#     x = 0:.1:10
#     n = 2
#     order = 2
#     lam =  -0.5
#     chi = 1.
#     psi = .15
#     a0 = 0
#     a = ones(n)
#     A = collect(I(n))
#     C = collect(I(n))
#     mu = zeros(n)
#     gam = zeros(n)
#
#     ccdf_true, pm_true = qfmgh(
#         x,
#         a0,
#         a,
#         A,
#         C,
#         mu,
#         gam,
#         lam,
#         chi,
#         psi;
#         do_spa=false
#         )
#
#     ccdf_spa, pm_spa = qfmgh(
#         x,
#         a0,
#         a,
#         A,
#         C,
#         mu,
#         gam,
#         lam,
#         chi,
#         psi;
#         do_spa=true,
#         order=order
#         )
#     plot(x, [ccdf_true, ccdf_spa])
# #end

#computes the pdf of the 2SLS estimator under multivariate t errors
using Random
Random.seed!(0)
order = 2
n = 25
beta = 0.
k = 1
Z = randn(n, k)
In = collect(I(n))
Z = In[:, 1:k] * (Z'*Z)^(.5) # distribution depends on Z only through Z'Z if R=I
alpha = 0.
if abs(alpha) >= 1.
     b = 0.
else
     b = 1 / sqrt(1-alpha^2)
end
R = Toeplitz(alpha .^(0:n-1), vcat(1., zeros(n-1))) * I
R[:, 1] = b * R[:, 1]
R = R ./ sqrt(tr(R * R') / n)
mu = .5
pi = [sqrt(mu / k)] * (Z' * Z)^(-.5) * ones(k, 1)
s2u = 1.
s2v = 1.
rho = 1.
suv = sqrt(s2u * s2v) * rho
nuvec = 3:2:15

h = 1e-6
Zp = Z * pi
Pz = Z * inv(Z' * Z) * Z'
xvec = beta-3:.02:beta+3
#E = eigen(.5*Pz+.5*Pz')
#[pp,ll]=eig(.5*Pz+.5*Pz')
nx = length(xvec)
pdf = zeros(length(nuvec), nx)
cdf = zeros(nx)
cdf2 = copy(cdf)

spapdf = zeros(length(nuvec), nx)
spacdf = zeros(nx)
spacdf2 = copy(cdf)

for nuloop=1:length(nuvec)
    nu = nuvec[nuloop]
    for loop = 1 : nx
        Si = ((nu-2) / nu * [s2u suv; suv s2v])^.5 # rescale to make Si^2 the covariance matrix
        S = kron(Si, R*R')
        b = xvec[loop] - beta
        a0 = (-b * (Zp' * Zp))[1]
        a = [Zp; -2*b*Zp]
        A = .5*[Pz*0 Pz; Pz -2*b*Pz]
        SAS = S * A * S
        SAS = .5*(SAS+SAS')
        E = eigen(SAS)
        P = E.vectors
        omega = E.values
        d = a' * S * P
        ccdf, _ = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=false)
        cdf[loop] = 1.0 - ccdf
        ccdf, _ = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=true, order=order)
        spacdf[loop] = 1.0 - ccdf
    end
    xvech = xvec .+ h
    for loop = 1 : nx
        Si = ((nu-2) / nu * [s2u suv; suv s2v])^.5 # rescale to make Si^2 the covariance matrix
        S = kron(Si, R*R')
        b = xvech[loop] - beta
        a0 = (-b * (Zp' * Zp))[1]
        a = [Zp; -2*b*Zp]
        A = .5*[Pz*0 Pz; Pz -2*b*Pz]
        SAS = S * A * S
        SAS = .5*(SAS+SAS')
        E = eigen(SAS)
        P = E.vectors
        omega = E.values
        d = a' * S * P
        ccdf, _ = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=false)
        cdf2[loop] = 1.0 .- ccdf
        ccdf, _ = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=true, order=order)
        spacdf2[loop] = 1.0 .- ccdf
    end
    pdf[nuloop, :] = (cdf2 .- cdf) ./ h
    spapdf[nuloop, :] = (spacdf2 .- spacdf) ./ h
end
spapdf[abs.(pdf.-spapdf) .> .1] .= NaN # filter errors due to nonsingularity around the mean
plot(xvec, pdf', color=collect(1:7)', labels=permutedims("v=".*string.(nuvec)))
plot!(xvec, spapdf', color=collect(1:7)', labels="", ls=:dash)
xlabel!("x")
ylabel!("pdf")
title!("Exact density (solid) and SPA (dashes)")
savefig("pdfs.svg")

module QuadraticFormsMGHyp

using LinearAlgebra
using Roots
using SpecialFunctions: besselk #, lgamma
using QuadGK: quadgk
using StatsFuns: normcdf, normpdf, logtwo
# work around https://github.com/JuliaMath/SpecialFunctions.jl/issues/186
# until https://github.com/JuliaDiff/ForwardDiff.jl/pull/419/ is merged
using Base.Math: libm
using ForwardDiff: Dual, value, partials, derivative
@inline lgamma(x::Float64) = ccall((:lgamma, libm), Float64, (Float64,), x)
@inline lgamma(x::Float32) = ccall((:lgammaf, libm), Float32, (Float32,), x)
@inline lgamma(d::Dual{T}) where {T} =
    Dual{T}(lgamma(value(d)), digamma(value(d)) * partials(d))
Base.@irrational rp 0.3183098861837906715 1 / big(pi)

export qfmgh

"""
    qfmgh(x, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)

Survivor function P(L>x) = 1-F(x) and
tail conditional mean, E(L|L>x), of

L = a0 + a' * X + X' * A * X, where:

   X = mu + W * gam + sqrt(W) * C * Z,
   Z ~ N(0, I), and W ~ GIG(lam, chi, psi), i.e.,
   X is distributed as multivariate GHyp.

Keyword arguments:

    `do_spa`: whether to return the exact result or a saddlepoint approximation
    `order`: order of the saddlepoint approximation

 (c) 2020 S.A. Broda
"""
function qfmgh end

qfmgh(x::Number, args...; kwargs...) = getindex.(qfmgh([x], args...; kwargs...), 1)

function qfmgh(
    x::AbstractVector,
    a0,
    a,
    A,
    C,
    mu,
    gam,
    lam,
    chi,
    psi;
    do_spa::Bool = false,
    order::Int = 2,
)
    if (do_spa && (lam>=0 || chi<= 0 ||  any(gam .!= 0 )))
        @warn "Saddlepoint approximation is inaccurate with these parameters."
    end

    CAC = Symmetric(C' * A * C)
    E = eigen(CAC)
    omega = E.values
    P = E.vectors
    muA = mu' * A
    CP = C * P
    gA = gam' * A
    c = a' * gam + 2 * muA * gam
    d = a' * CP + 2 * muA * CP
    e = 2 * gA * CP
    de = d .* e
    d2 = d .* d
    e2 = e .* e
    k = gA * gam
    kk = a0 + a' * mu + muA * mu
    LK2 = lklam(lam, chi, psi)
    qq = x .- kk
    ccdf = similar(float(x))
    pm = similar(float(x))
    M, alpha2p, dM0da1, alpha1p, M0, lrhop = get_funcs(omega, de, e2, d2, c, k, LK2, lam, chi, psi)
    #Threads.@threads
    for i = 1:length(qq)
        q = qq[i]
        if ~do_spa
            ccdf[i], _ = quadgk(s -> imag(M(1im * s, -1im * s * q) / s), 0.0, Inf)
            ccdf[i] = M(0, 0) / 2 + rp * ccdf[i]
            M2(s, t) = M(s, -q * s) * alpha2p(s) + dM0da1(s, -q * s) * alpha1p(s) + M0(s, -q * s) * lrhop(s)
            pm[i], _ = quadgk(s -> imag(M2(1im * s, -1im * s * q) / s), 0.0, Inf)
            pm[i] = (M2(0, 0) / 2 + rp * pm[i]) / ccdf[i] + kk
        else
            ccdf[i] = compute_spa(s -> 1, s -> log(M(s, -q * s)), order)
            pm[i] =
                (
                    compute_spa(alpha2p, s -> log(M(s, -q * s)), order) +
                    compute_spa(alpha1p, s -> log(dM0da1(s, -q * s)), order) +
                    compute_spa(lrhop, s -> log(M0(s, -q * s)), order)
                ) / ccdf[i] + kk
        end
    end
    return ccdf, pm
end

@inline function lklam(lam, chi, psi)
    if chi == 0
        if real(psi) < 0
            return Inf
        else
            return -lam * log(psi / 2) + lgamma(lam)
        end
    elseif psi == 0
        if real(chi) < 0
            return Inf
        else
            return lam * log(chi / 2) + lgamma(-lam)
        end
    elseif real(chi * psi) < 0
        return Inf
    else
        scp = sqrt(chi * psi)
        return logtwo + lam * log(chi / psi) / 2 + log(besselk(lam, scp))
    end
end

function compute_spa(g, h, order)
    h0 = h(0.0)
    g0 = g(0.0)
    hp = s -> derivative(h, s)
    hpp = s -> derivative(hp, s)
    shat = find_zero(hp, 0.0)
    what = sign(shat) * sqrt(-2 * (h(shat) - h0))
    H = hpp(shat)
    uhat = shat * sqrt(H)
    ghat = g(shat)
    spa = exp(h0) * (g0 * (1 - normcdf(what)) + normpdf(what) * (ghat / uhat - g0 / what))
    if order == 2
        gp = s -> derivative(g, s)
        gpp = s -> derivative(gp, s)
        hppp = s -> derivative(hpp, s)
        hpppp = s -> derivative(hppp, s)
        k3 = hppp(shat) / H^(3 // 2)
        k4 = hpppp(shat) / H^2
        T1 = ghat / uhat * ((k4 / 8 - 5 * k3^2 / 24) - uhat^-2 - k3 / (2 * uhat))
        T2 = shat * gp(shat) / uhat * (1 / uhat^2 + k3 / (2 * uhat))
        T3 = -gpp(shat) * shat^2 / (2 * uhat^3) + g0 * what^-3
        spa = spa + exp(h0) * normpdf(what) * (T1 + T2 + T3)
    end
    return spa
end

function get_funcs(omega, de, e2, d2, c, k, LK2, lam, chi, psi)
    M(s, t) = exp(
        lklam(
            lam,
            chi - 2 * (0.5 * s^2 * sum(d2 ./ (1 .- 2 * omega * s)) + t),
            psi - 2 * (k * s + 0.5 * s^2 * sum(e2 ./ (1 .- 2 * omega * s))),
        ) - LK2 +
        s * c +
        s^2 * sum(de ./ (1 .- 2 * omega * s)) +
        0.5 * sum(log, 1 ./ (1 .- 2 * omega * s)),
    )
    M0(s, t) = exp(
        lklam(
            lam + 1,
            chi - 2 * (0.5 * s^2 * sum(d2 ./ (1 .- 2 * omega * s)) + t),
            psi - 2 * (k * s + 0.5 * s^2 * sum(e2 ./ (1 .- 2 * omega * s))),
        ) - LK2 +
        s * c +
        s^2 * sum(de ./ (1 .- 2 * omega * s)) +
        0.5 * sum(log, 1 ./ (1 .- 2 * omega * s)),
    )
    dM0da1(s, t) = exp(
        lklam(
            lam + 1,
            chi - 2 * (0.5 * s^2 * sum(d2 ./ (1 .- 2 * omega * s)) + t),
            psi - 2 * (k * s + 0.5 * s^2 * sum(e2 ./ (1 .- 2 * omega * s))),
        ) - LK2 +
        s * c +
        s^2 * sum(de ./ (1 .- 2 * omega * s)) +
        0.5 * sum(log, 1 ./ (1 .- 2 * omega * s)),
    )
    lrhop(s) =
        c +
        2 * s * sum(de ./ (1 .- 2 * omega * s)) +
        2 * s^2 * sum(de .* omega ./ (1 .- 2 * omega * s) .^ 2) +
        sum(omega ./ (1 .- 2 * omega * s))
    alpha1p(s) =
        k +
        s * sum(e2 ./ (1 .- 2 * omega * s)) +
        s^2 * sum(e2 .* omega ./ (1 .- 2 * omega * s) .^ 2)
    alpha2p(s) =
        s * sum(d2 ./ (1 .- 2 * omega * s)) +
        s^2 * sum(d2 .* omega ./ (1 .- 2 * omega * s) .^ 2)
    return M, alpha2p, dM0da1, alpha1p, M0, lrhop
end
end

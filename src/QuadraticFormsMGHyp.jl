module QuadraticFormsMGHyp

using LinearAlgebra
using Roots
using SpecialFunctions: besselk #, lgamma
using QuadGK: quadgk
using StatsFuns: normcdf, normpdf, logtwo
using PrecompileTools
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
    qfmgh(x::Union{AbstractVector{<:Real}, Real}, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)

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

qfmgh(x::Real, args...; kwargs...) = getindex.(qfmgh([x], args...; kwargs...), 1)

function qfmgh(
    x::AbstractVector{<:Real},
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
    lM, alpha2p, ldM0da1, alpha1p, lM0, lrhop = get_funcs(omega, de[:], e2[:], d2[:], c, k, LK2, lam, chi, psi)
    shat0 = zeros(Threads.nthreads())
    shat2 = zeros(Threads.nthreads())
    shat3 = zeros(Threads.nthreads())
    Threads.@threads for i = 1:length(qq)
        q = qq[i]
        if ~do_spa
            ccdf[i], _ = quadgk(s -> imag(exp(lM(1im * s, -1im * s * q)) / s), 0.0, Inf)
            ccdf[i] = exp(lM(0, 0)) / 2 + rp * ccdf[i]
            M2(s, t) = exp(lM(s, -q * s)) * alpha2p(s) + exp(ldM0da1(s, -q * s)) * alpha1p(s) + exp(lM0(s, -q * s)) * lrhop(s)
            pm[i], _ = quadgk(s -> imag(M2(1im * s, -1im * s * q) / s), 0.0, Inf)
            pm[i] = (M2(0, 0) / 2 + rp * pm[i]) / ccdf[i] + kk
        else
            ccdf[i], shat0[Threads.threadid()] = compute_spa(s -> 1, s -> lM(s, -q * s), order, shat0[Threads.threadid()])
            I1 = all(d.==0) ? 0. : compute_spa(alpha2p, s -> lM(s, -q * s), order, shat0[Threads.threadid()], false)[1]
            I2, shat2[Threads.threadid()] = all(gam.==0) ? (0., shat2[Threads.threadid()]) : compute_spa(alpha1p, s -> ldM0da1(s, -q * s), order, shat2[Threads.threadid()])
            I3, shat3[Threads.threadid()] = compute_spa(lrhop, s -> lM0(s, -q * s), order, shat3[Threads.threadid()])
            pm[i] = (I1 + I2 + I3) / ccdf[i] + kk
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

function compute_spa(g, h, order, shat=0., solve=true)
    h0 = h(0.0)
    g0 = g(0.0)
    hp = s -> derivative(h, s)
    hpp = s -> derivative(hp, s)

    if solve # otherwise, keep starting value
        shat = find_zero(hp, shat)
    end

    if abs(h(shat) - h0) < 1e-5
        @warn("Saddlepoint approximation is inaccurate; returning NaN.")
        spa, shat = NaN, 0.
    else
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
    end
    return spa, shat
end

function get_funcs(omega, de, e2, d2, c, k, LK2, lam, chi, psi)
    @inline logTheta(s, t, j) = @fastmath (t1 = 0.; t2 = 0.; t3 = 0.; t4 = 0.;
                              @inbounds @simd ivdep for i = 1:length(omega)
                                  nu = 1 / (1 - 2 * omega[i] * s)
                                  t1 += d2[i] * nu
                                  t2 += e2[i] * nu
                                  t3 += de[i] * nu
                                  t4 += log(nu)
                              end;
                              lklam(lam + j,
                                        chi - 2 * (0.5 * s^2 * t1 + t),
                                        psi - 2 * (k * s + 0.5 * s^2 * t2),
                                        ) - LK2 + s * c + s^2 * t3 + 0.5 * t4,
                                    )
    lM(s, t) = logTheta(s, t, 0)
    lM0(s, t) = logTheta(s, t, 1)
    ldM0da1(s, t) = logTheta(s, t, 2)
    lrhop(s) = @fastmath (t=0.; @inbounds @simd ivdep for i = 1:length(omega)
                        nu = 1 / (1 - 2 * omega[i] * s)
                        t += 2 * s * de[i] * nu + 2 * s^2 * de[i] * omega[i] * nu +  omega[i] * nu
                      end;
                      t + c
                )
    alpha1p(s) = @fastmath (t=0.; @inbounds @simd ivdep for i = 1:length(omega)
                          nu = 1 / (1 - 2 * omega[i] * s)
                          t += s * e2[i] * nu + s^2 * e2[i] * omega[i] * nu^2
                        end;
                        t + k
                )

    alpha2p(s) = @fastmath (t=0.; @inbounds @simd ivdep for i = 1:length(omega)
                          nu = 1 / (1 - 2 * omega[i] * s)
                          t += s * d2[i] * nu + s^2 * d2[i] * omega[i] * nu^2
                        end;
                        t
                 )
    return lM, alpha2p, ldM0da1, alpha1p, lM0, lrhop
end

@static if VERSION >= v"1.9.0-alpha1"
	@compile_workload begin
		ccdf, pm = qfmgh(5.991, 0., zeros(2), [1.0 0.; 0. 1.], [1. 0.; 0. 1.], zeros(2), zeros(2), -10/2, 10, 0., do_spa=false)
        ccdf, pm = qfmgh(5.991, 0., zeros(2), [1.0 0.; 0. 1.], [1. 0.; 0. 1.], zeros(2), zeros(2), -10/2, 10, 0., do_spa=true)
	end # precompile block
end # if
end # module

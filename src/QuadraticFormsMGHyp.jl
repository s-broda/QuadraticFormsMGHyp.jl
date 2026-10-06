module QuadraticFormsMGHyp

using LinearAlgebra
using Roots
using SpecialFunctions: besselk, besselkx #, lgamma
using StatsFuns: normcdf, normpdf, logtwo
using PrecompileTools: @compile_workload
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

The exact result is a mapped Gauss–Legendre rule. The order doubles until successive refinements agree to a relative tolerance of 1e-7, up to 4096 nodes. Lists of more than 24 thresholds are reduced to a Chebyshev series. `chi = psi = Inf` is the Gaussian limit: the mixer is the constant 1, for any `lam`. A low-rank Gaussian spectrum is integrated in panels, because its characteristic function decays only as a power of the frequency.

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
    gaussian = isinf(chi) && isinf(psi) && chi > 0 && psi > 0
    if (do_spa && (lam>=0 || chi<= 0 ||  any(gam .!= 0 ) || gaussian))
        @warn "Saddlepoint approximation is inaccurate with these parameters."
    end
    if !do_spa || gaussian
        return quadrature_eval(x, a0, a, A, C, mu, gam, lam, chi, psi)
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
    Threads.@threads for i = 1:length(qq)
        if !haskey(task_local_storage(), :shat0)
            task_local_storage(:shat0, 0.)
            task_local_storage(:shat2, 0.)
            task_local_storage(:shat3, 0.)
        end
        q = qq[i]
        ccdf[i], shat0 = compute_spa(s -> 1, s -> lM(s, -q * s), order, task_local_storage(:shat0))
            task_local_storage(:shat0, shat0)
            I1 = all(d.==0) ? 0. : compute_spa(alpha2p, s -> lM(s, -q * s), order, task_local_storage(:shat0), false)[1]
            I2, shat2 = all(gam.==0) ? (0., task_local_storage(:shat2)) : compute_spa(alpha1p, s -> ldM0da1(s, -q * s), order, task_local_storage(:shat2))
            task_local_storage(:shat2, shat2)
            I3, shat3 = compute_spa(lrhop, s -> lM0(s, -q * s), order, task_local_storage(:shat3))
            task_local_storage(:shat3, shat3)
            pm[i] = (I1 + I2 + I3) / ccdf[i] + kk
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
                                        ) - LK2 + s * c + s^2 * t3 + 0.5 * t4
                                    )
    lM(s, t) = logTheta(s, t, 0)
    lM0(s, t) = logTheta(s, t, 1)
    ldM0da1(s, t) = logTheta(s, t, 2)
    lrhop(s) = @fastmath (t=0.; @inbounds @simd ivdep for i = 1:length(omega)
                        nu = 1 / (1 - 2 * omega[i] * s)
                        # ν^2 on the s^2 term, ν on the last term: (log ρ)' / i.
                        t += 2 * s * de[i] * nu + 2 * s^2 * de[i] * omega[i] * nu * nu + omega[i] * nu
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

include("quadrature.jl")

@static if VERSION >= v"1.9.0-alpha1"
    # One cheap call per compiled body, inlined so each specialization is a
    # precompile root. A `let` or a non-constant global hides that root.
    # Tails that climb toward 4096 nodes use the same methods and are omitted.
    @compile_workload begin
        # NIG, direct and Chebyshev. A is Diagonal, as in the documented example.
        qfmgh(5.0, 0.0, zeros(10), 0.5I(10), Matrix{Float64}(I, 10, 10), zeros(10), zeros(10), -0.5, 1.0, 1.0)
        qfmgh(range(3.5, 17.5; length=32), 0.0, zeros(10), 0.5I(10), Matrix{Float64}(I, 10, 10), zeros(10), zeros(10), -0.5, 1.0, 1.0)
        # ψ = 0 closed form, half-integer Bessel, generic Bessel.
        qfmgh(5.991, 0.0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -5.0, 10.0, 0.0)
        qfmgh(1.0, 0.0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -1.5, 2.0, 1.5)
        qfmgh(1.0, 0.0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -1.3, 2.0, 1.1)
        # Gaussian: mapped linear, low-rank panels, and the constant atom.
        qfmgh(1.0, 0.0, [1.0], zeros(1, 1), ones(1, 1), zeros(1), zeros(1), 0.0, Inf, Inf)
        qfmgh(10.0, 0.0, zeros(1), ones(1, 1), ones(1, 1), zeros(1), zeros(1), 0.0, Inf, Inf)
        qfmgh(0.0, 0.0, zeros(1), ones(1, 1), ones(1, 1), zeros(1), zeros(1), 0.0, Inf, Inf)
        qfmgh([-1.0, 1.0], 0.0, zeros(1), zeros(1, 1), ones(1, 1), zeros(1), zeros(1), 0.0, Inf, Inf)
        # Saddlepoint. These parameters do not trip the inaccuracy warning.
        qfmgh(5.991, 0.0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -5.0, 10.0, 0.0; do_spa=true)
        # Integer scalars. `qfmgh(1, 0, …, -1, 2, 1)` and an integer threshold
        # with otherwise floating arguments are separate specializations.
        qfmgh(1, 0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -1, 2, 1)
        qfmgh(1, 0.0, zeros(2), [1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], zeros(2), zeros(2), -1.5, 2.0, 1.5)
    end
end # if
end # module

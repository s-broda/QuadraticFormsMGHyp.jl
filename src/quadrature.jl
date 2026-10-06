# Mapped Gauss–Legendre evaluation of the Gil-Pelaez integral.
# The node count starts at the law's default and doubles until successive
# refinements agree to relative tolerance 1e-7, or until 4096 nodes.
# The truncation and the Chebyshev growth match es4mgh.c.

const KIND_GENERAL = 0
const KIND_NIG = 1
const KIND_HALF = 2
const KIND_PSI0 = 3
const KIND_GAUSS = 4
const IM = ComplexF64(0.0, 1.0)
const LOG2F = log(2.0)
const REL_QUAD = 1e-7
const NNODE_CAP = 4096
# A Gaussian spectrum whose mapped truncation sits above this bound decays
# too slowly for one panel. The panel rule below is used instead.
const GAUSS_FAST_UB = 0.91
const GAUSS_XT = 800.0

struct Prepared
    kind::Int
    kk::Float64
    lam::Float64
    LK2::Float64
    M20::Float64
    need_a1::Bool
    nnode::Int
    u::Vector{Float64}
    wo::Vector{Float64}
    chi_base::Vector{ComplexF64}
    psi_node::Vector{ComplexF64}
    lrho::Vector{ComplexF64}
    a2p::Vector{ComplexF64}
    a1p::Vector{ComplexF64}
    lrp::Vector{ComplexF64}
    log_psi::Vector{ComplexF64}
end

# Spectral reduction, reused while the node count is doubled.
# nlo is the coarse order that agreed with its refinement.
mutable struct Layout
    kind::Int
    kk::Float64
    lam::Float64
    chi::Float64
    psi::Float64
    LK2::Float64
    M20::Float64
    need_a1::Bool
    ccoef::Float64
    k::Float64
    D2z::Float64
    E2z::Float64
    DEz::Float64
    ub::Float64
    omega::Vector{Float64}
    d2::Vector{Float64}
    e2::Vector{Float64}
    de::Vector{Float64}
    nlo::Int
end

function legendre_pd(n::Int, x::Float64)
    n == 0 && return 1.0, 0.0
    n == 1 && return x, 1.0
    p0 = 1.0
    p1 = x
    @inbounds for k = 2:n
        p = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k
        p0 = p1
        p1 = p
    end
    d = n * (x * p1 - p0) / (x * x - 1.0)
    return p1, d
end

function gauss_legendre(n::Int)
    x = Vector{Float64}(undef, n)
    w = Vector{Float64}(undef, n)
    @inbounds for i = 1:n
        theta = π * (i - 1 + 0.75) / (n + 0.5)
        xi = cos(theta)
        for _ = 1:20
            Pn, dPn = legendre_pd(n, xi)
            dx = Pn / dPn
            xi -= dx
            abs(dx) < 1e-15 && break
        end
        _, dPn = legendre_pd(n, xi)
        x[i] = xi
        w[i] = 2.0 / ((1.0 - xi * xi) * dPn * dPn)
    end
    return x, w
end

const GAUSS_GX, GAUSS_GW = gauss_legendre(12)

function integration_ub(omega::AbstractVector{Float64})
    o = sort!(abs.(filter(a -> a > 0.0, omega)), rev = true)
    m = length(o)
    epsabs = 1e-10
    ub = 1.0
    logsum = 0.0
    @inbounds for i = 1:m
        logsum += log(2.0 * o[i])
        ubnew = (-1.0 / i) * (2.0 * log(π) + 2.0 * log(i) + 2.0 * log(epsabs) - log(4.0) + logsum)
        ubnew = exp(ubnew)
        ubnew = sqrt(ubnew) / (1.0 + sqrt(ubnew))
        ubnew < ub && (ub = ubnew)
    end
    !(ub > 1e-6 && ub < 1.0) && (ub = 0.99)
    ub > 0.999999 && (ub = 0.999999)
    return ub
end

function half_order(nu::Float64)
    a = abs(nu)
    h = round(a - 0.5)
    return h >= -0.5 && h < 60.0 && abs(a - (h + 0.5)) < 1e-10
end

function default_nnode(kind::Int)
    (kind == KIND_NIG || kind == KIND_GENERAL || kind == KIND_GAUSS) && return 32
    kind == KIND_HALF && return 48
    return 64
end

function cheb_order(kind::Int)
    # First degree tried. The normal-inverse-Gaussian tail has settled by 20.
    # Other laws start higher and still grow until the last coefficient is small.
    (kind == KIND_NIG || kind == KIND_GAUSS) && return 20
    return 48
end

function lklam_real_f(lam::Float64, chi::Float64, psi::Float64)
    chi == 0.0 && return -lam * log(psi * 0.5) + lgamma(lam)
    psi == 0.0 && return lam * log(chi * 0.5) + lgamma(-lam)
    z = sqrt(chi * psi)
    # besselkx is exp(z) K_λ(z). Subtract z to recover log K, as log_scaled_K does.
    logK = log(besselkx(lam, z)) - z
    return LOG2F + 0.5 * lam * log(chi / psi) + logK
end

function lklam_c(lam::Float64, chi::ComplexF64, psi::ComplexF64)
    if chi == 0
        real(psi) < 0 && return ComplexF64(Inf)
        return -lam * log(psi * 0.5) + lgamma(lam)
    end
    if psi == 0
        real(chi) < 0 && return ComplexF64(Inf)
        return lam * log(chi * 0.5) + lgamma(-lam)
    end
    prod = chi * psi
    if real(prod) < 0 && abs(imag(prod)) <= 1e-14 * (1.0 + abs(real(prod)))
        return ComplexF64(Inf)
    end
    z = sqrt(prod)
    logK = log(besselkx(lam, z)) - z
    return LOG2F + 0.5 * lam * log(chi / psi) + logK
end

@inline function imag_exp(phase::ComplexF64)
    re = real(phase)
    (re < -700.0 || re > 700.0) && return 0.0
    return exp(re) * sin(imag(phase))
end

@inline function imag_exp_mul(phase::ComplexF64, c::ComplexF64)
    re = real(phase)
    (re < -700.0 || re > 700.0) && return 0.0
    e = exp(re)
    s = sin(imag(phase))
    co = cos(imag(phase))
    return e * (s * real(c) + co * imag(c))
end

function prepare(a0, a, A, C, mu, gam, lam, chi, psi)
    a = collect(Float64, vec(a))
    d = length(a)
    A = reshape(collect(Float64, vec(A)), d, d)
    C = reshape(collect(Float64, vec(C)), d, d)
    mu = collect(Float64, vec(mu))
    gam = collect(Float64, vec(gam))
    lam = Float64(lam)
    chi = Float64(chi)
    psi = Float64(psi)
    As = 0.5 .* (A .+ transpose(A))
    F = eigen(Symmetric(transpose(C) * As * C))
    omega = F.values
    P = F.vectors
    muA = As * mu
    gA = As * gam
    ccoef = dot(a, gam) + 2.0 * dot(muA, gam)
    kk = Float64(a0) + dot(a, mu) + dot(muA, mu)
    k = dot(gA, gam)
    CP = C * P
    dvec = transpose(CP) * (a .+ 2.0 .* muA)
    evec = transpose(CP) * (2.0 .* gA)
    return prepare_spectral(omega, dvec, evec, ccoef, k, kk, lam, chi, psi)
end

function prepare_spectral(omega_all, d_all, e_all, ccoef, k, kk, lam, chi, psi)
    omega_all = collect(Float64, vec(omega_all))
    d_all = collect(Float64, vec(d_all))
    e_all = collect(Float64, vec(e_all))
    ne_all = length(omega_all)
    omax = isempty(omega_all) ? 0.0 : maximum(abs, omega_all)
    otol = 1e-12 * (omax > 1.0 ? omax : 1.0)
    omega = Float64[]
    d2 = Float64[]
    e2 = Float64[]
    de = Float64[]
    D2z = 0.0
    E2z = 0.0
    DEz = 0.0
    @inbounds for i = 1:ne_all
        di = d_all[i]
        ei = e_all[i]
        if abs(omega_all[i]) <= otol
            D2z += di * di
            E2z += ei * ei
            DEz += di * ei
        else
            push!(omega, omega_all[i])
            push!(d2, di * di)
            push!(e2, ei * ei)
            push!(de, di * ei)
        end
    end
    ne = length(omega)
    need_a1 = abs(k) > 0.0 || E2z > 0.0
    @inbounds for i = 1:ne
        (e2[i] != 0.0 || de[i] != 0.0) && (need_a1 = true)
    end
    ub = integration_ub(omega_all)
    sum_om = isempty(omega_all) ? 0.0 : sum(omega_all)
    # χ = ψ = +∞ is the degenerate mixer W ≡ 1, for any λ.
    # X is then Gaussian with mean μ + γ and covariance CCᵀ.
    # The truncation uses |ω|: a negative eigenvalue decays like a positive one.
    if isinf(chi) && isinf(psi) && chi > 0.0 && psi > 0.0
        ub = integration_ub(abs.(omega_all))
        M20 = k + ccoef + sum_om
        return Layout(KIND_GAUSS, Float64(kk), Float64(lam), Float64(chi), Float64(psi),
                      0.0, M20, need_a1, Float64(ccoef), Float64(k),
                      D2z, E2z, DEz, ub, omega, d2, e2, de, 0)
    end
    LK2 = lklam_real_f(lam, chi, psi)
    kind = KIND_GENERAL
    if abs(lam + 0.5) < 1e-12 && chi > 0.0 && psi > 0.0
        kind = KIND_NIG
    elseif psi == 0.0 && abs(k) == 0.0 && E2z == 0.0
        kind = KIND_PSI0
    elseif half_order(lam)
        kind = KIND_HALF
    end
    if kind == KIND_PSI0 && need_a1
        kind = half_order(lam) ? KIND_HALF : KIND_GENERAL
    end
    lm1 = lklam_real_f(lam + 1.0, chi, psi) - LK2
    lm2 = lklam_real_f(lam + 2.0, chi, psi) - LK2
    # The skewness term carries E[W^2]. When that coefficient is zero the
    # product is zero even if the moment is infinite (integer λ, ψ = 0).
    skew = k == 0.0 ? 0.0 : real(exp(lm2) * k)
    M20 = skew + real(exp(lm1) * (ccoef + sum_om))
    return Layout(kind, Float64(kk), Float64(lam), Float64(chi), Float64(psi),
                  LK2, M20, need_a1, Float64(ccoef), Float64(k),
                  D2z, E2z, DEz, ub, omega, d2, e2, de, 0)
end

function materialize(L::Layout, nn::Int)
    gx, gw = gauss_legendre(nn)
    u = Vector{Float64}(undef, nn)
    wo = Vector{Float64}(undef, nn)
    chi_base = Vector{ComplexF64}(undef, nn)
    psi_node = Vector{ComplexF64}(undef, nn)
    lrho = Vector{ComplexF64}(undef, nn)
    a2pv = Vector{ComplexF64}(undef, nn)
    a1pv = Vector{ComplexF64}(undef, nn)
    lrp = Vector{ComplexF64}(undef, nn)
    log_psi = Vector{ComplexF64}(undef, nn)
    ne = length(L.omega)
    ub = L.ub
    @inbounds for i = 1:nn
        v = 0.5 * ub * (gx[i] + 1.0)
        wv = 0.5 * ub * gw[i]
        ss = v / (1.0 - v)
        uu = ss * ss
        dudv = 2.0 * ss / ((1.0 - v) * (1.0 - v))
        u[i] = uu
        wo[i] = wv * dudv / uu
        s = IM * uu
        s2 = s * s
        t1 = zero(ComplexF64)
        t2 = zero(ComplexF64)
        t3 = zero(ComplexF64)
        t4 = zero(ComplexF64)
        a2p = zero(ComplexF64)
        a1p = zero(ComplexF64)
        lr = zero(ComplexF64)
        for j = 1:ne
            nu = 1.0 / (1.0 - 2.0 * L.omega[j] * s)
            nu2 = nu * nu
            t1 += L.d2[j] * nu
            t2 += L.e2[j] * nu
            t3 += L.de[j] * nu
            t4 += log(nu)
            a2p += s * L.d2[j] * nu + s2 * L.d2[j] * L.omega[j] * nu2
            a1p += s * L.e2[j] * nu + s2 * L.e2[j] * L.omega[j] * nu2
            lr += 2.0 * s * L.de[j] * nu + 2.0 * s2 * L.de[j] * L.omega[j] * nu2 + L.omega[j] * nu
        end
        t1 += L.D2z
        t2 += L.E2z
        t3 += L.DEz
        a2p += s * L.D2z
        a1p += s * L.E2z + L.k
        lr += 2.0 * s * L.DEz + L.ccoef
        if L.kind == KIND_GAUSS
            # log ρ + α₁ + α₂. Evaluation subtracts i u q, which is i t.
            chi_base[i] = 0.0
            psi_node[i] = 0.0
            lrho[i] = s * L.ccoef + s2 * t3 + 0.5 * t4 + L.k * s + 0.5 * s2 * (t1 + t2)
            a2pv[i] = a2p
            a1pv[i] = a1p
            lrp[i] = lr
            log_psi[i] = 0.0
            continue
        end
        chi_base[i] = L.chi - s2 * t1
        pnode = L.psi - 2.0 * (L.k * s + 0.5 * s2 * t2)
        psi_node[i] = pnode
        lrho[i] = s * L.ccoef + s2 * t3 + 0.5 * t4
        a2pv[i] = a2p
        a1pv[i] = a1p
        lrp[i] = lr
        pn = real(pnode)^2 + imag(pnode)^2
        if pn == 0.0
            log_psi[i] = 0.0
        else
            log_psi[i] = log(pnode)
        end
    end
    return Prepared(L.kind, L.kk, L.lam, L.LK2, L.M20, L.need_a1, nn,
                    u, wo, chi_base, psi_node, lrho, a2pv, a1pv, lrp, log_psi)
end

function finish_tail(E::Prepared, Ic::Float64, Ip::Float64)
    cval = 0.5 + rp * Ic
    es = (0.5 * E.M20 + rp * Ip) / cval + E.kk
    return cval, es
end

function eval_gauss!(E::Prepared, x, ccdf, es)
    nn = E.nnode
    need = E.need_a1
    @inbounds for qi in eachindex(x)
        q = x[qi] - E.kk
        Ic = 0.0
        Ip = 0.0
        for i = 1:nn
            lm = E.lrho[i] - IM * (E.u[i] * q)
            Ic += E.wo[i] * imag_exp(lm)
            acc = imag_exp_mul(lm, E.a2p[i]) + imag_exp_mul(lm, E.lrp[i])
            if need
                acc += imag_exp_mul(lm, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
        ccdf[qi], es[qi] = finish_tail(E, Ic, Ip)
    end
    return nothing
end

function eval_nig!(E::Prepared, x, ccdf, es)
    nn = E.nnode
    LK2 = E.LK2
    log_half_pi = log(0.5 * π)
    need = E.need_a1
    @inbounds for qi in eachindex(x)
        q = x[qi] - E.kk
        twoq = 2.0 * q
        Ic = 0.0
        Ip = 0.0
        for i = 1:nn
            chi = E.chi_base[i] + IM * (E.u[i] * twoq)
            z = sqrt(chi * E.psi_node[i])
            logK = -z + 0.5 * (log_half_pi - log(z))
            lrat = log(chi) - E.log_psi[i]
            lm0 = LOG2F - 0.25 * lrat + logK - LK2 + E.lrho[i]
            Ic += E.wo[i] * imag_exp(lm0)
            lm1 = LOG2F + 0.25 * lrat + logK - LK2 + E.lrho[i]
            acc = imag_exp_mul(lm0, E.a2p[i]) + imag_exp_mul(lm1, E.lrp[i])
            if need
                lm2 = LOG2F + 0.75 * lrat + logK + log(1.0 + 1.0 / z) - LK2 + E.lrho[i]
                acc += imag_exp_mul(lm2, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
        ccdf[qi], es[qi] = finish_tail(E, Ic, Ip)
    end
    return nothing
end

function eval_psi!(E::Prepared, x, ccdf, es)
    lam = E.lam
    LK2 = E.LK2
    c0 = lgamma(-lam) - lam * LOG2F - LK2
    c1 = lgamma(-(lam + 1.0)) - (lam + 1.0) * LOG2F - LK2
    c2 = lgamma(-(lam + 2.0)) - (lam + 2.0) * LOG2F - LK2
    nn = E.nnode
    need = E.need_a1
    @inbounds for qi in eachindex(x)
        q = x[qi] - E.kk
        Ic = 0.0
        Ip = 0.0
        for i = 1:nn
            chi = E.chi_base[i] + IM * (2.0 * E.u[i] * q)
            lc = log(chi)
            base = E.lrho[i]
            lm0 = lam * lc + c0 + base
            Ic += E.wo[i] * imag_exp(lm0)
            lm1 = (lam + 1.0) * lc + c1 + base
            acc = imag_exp_mul(lm0, E.a2p[i]) + imag_exp_mul(lm1, E.lrp[i])
            if need
                lm2 = (lam + 2.0) * lc + c2 + base
                acc += imag_exp_mul(lm2, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
        ccdf[qi], es[qi] = finish_tail(E, Ic, Ip)
    end
    return nothing
end

function eval_point(E::Prepared, q::Float64)
    nn = E.nnode
    lam = E.lam
    LK2 = E.LK2
    Ic = 0.0
    Ip = 0.0
    need = E.need_a1
    if E.kind == KIND_GAUSS
        @inbounds for i = 1:nn
            lm = E.lrho[i] - IM * (E.u[i] * q)
            Ic += E.wo[i] * imag_exp(lm)
            acc = imag_exp_mul(lm, E.a2p[i]) + imag_exp_mul(lm, E.lrp[i])
            if need
                acc += imag_exp_mul(lm, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
    elseif E.kind == KIND_NIG
        log_half_pi = log(0.5 * π)
        @inbounds for i = 1:nn
            chi = E.chi_base[i] + IM * (2.0 * E.u[i] * q)
            z = sqrt(chi * E.psi_node[i])
            logK = -z + 0.5 * (log_half_pi - log(z))
            lrat = log(chi) - E.log_psi[i]
            lm0 = LOG2F - 0.25 * lrat + logK - LK2 + E.lrho[i]
            Ic += E.wo[i] * imag_exp(lm0)
            lm1 = LOG2F + 0.25 * lrat + logK - LK2 + E.lrho[i]
            acc = imag_exp_mul(lm0, E.a2p[i]) + imag_exp_mul(lm1, E.lrp[i])
            if need
                lm2 = LOG2F + 0.75 * lrat + logK + log(1.0 + 1.0 / z) - LK2 + E.lrho[i]
                acc += imag_exp_mul(lm2, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
    elseif E.kind == KIND_PSI0
        c0 = lgamma(-lam) - lam * LOG2F - LK2
        c1 = lgamma(-(lam + 1.0)) - (lam + 1.0) * LOG2F - LK2
        c2 = lgamma(-(lam + 2.0)) - (lam + 2.0) * LOG2F - LK2
        @inbounds for i = 1:nn
            chi = E.chi_base[i] + IM * (2.0 * E.u[i] * q)
            lc = log(chi)
            base = E.lrho[i]
            lm0 = lam * lc + c0 + base
            Ic += E.wo[i] * imag_exp(lm0)
            lm1 = (lam + 1.0) * lc + c1 + base
            acc = imag_exp_mul(lm0, E.a2p[i]) + imag_exp_mul(lm1, E.lrp[i])
            if need
                lm2 = (lam + 2.0) * lc + c2 + base
                acc += imag_exp_mul(lm2, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
    else
        @inbounds for i = 1:nn
            chi = E.chi_base[i] + IM * (2.0 * E.u[i] * q)
            psi = E.psi_node[i]
            base = -LK2 + E.lrho[i]
            lm0 = lklam_c(lam, chi, psi) + base
            Ic += E.wo[i] * imag_exp(lm0)
            lm1 = lklam_c(lam + 1.0, chi, psi) + base
            acc = imag_exp_mul(lm0, E.a2p[i]) + imag_exp_mul(lm1, E.lrp[i])
            if need
                lm2 = lklam_c(lam + 2.0, chi, psi) + base
                acc += imag_exp_mul(lm2, E.a1p[i])
            end
            Ip += E.wo[i] * acc
        end
    end
    return finish_tail(E, Ic, Ip)
end

function eval_direct!(E::Prepared, x, ccdf, es)
    if E.kind == KIND_GAUSS
        eval_gauss!(E, x, ccdf, es)
        return nothing
    elseif E.kind == KIND_NIG
        eval_nig!(E, x, ccdf, es)
        return nothing
    elseif E.kind == KIND_PSI0
        eval_psi!(E, x, ccdf, es)
        return nothing
    end
    n = length(x)
    if n >= 4 && Threads.nthreads() > 1
        Threads.@threads for i = 1:n
            @inbounds ccdf[i], es[i] = eval_point(E, x[i] - E.kk)
        end
    else
        @inbounds for i = 1:n
            ccdf[i], es[i] = eval_point(E, x[i] - E.kk)
        end
    end
    return nothing
end

function cheb_coeffs!(a, yq)
    m = length(yq)
    @inbounds for k = 0:m-1
        s = 0.0
        for j = 0:m-1
            s += yq[j + 1] * cos(π * k * (j + 0.5) / m)
        end
        a[k + 1] = k == 0 ? s / m : 2.0 * s / m
    end
    return a
end

function clenshaw(a, t::Float64)
    m = length(a)
    u2 = 0.0
    u1 = 0.0
    @inbounds for k = m:-1:2
        u = 2.0 * t * u1 - u2 + a[k]
        u2 = u1
        u1 = u
    end
    return t * u1 - u2 + a[1]
end

function interp!(y, x, yq, mid::Float64, half::Float64)
    a = Vector{Float64}(undef, length(yq))
    cheb_coeffs!(a, yq)
    invh = half == 0.0 ? 0.0 : 1.0 / half
    @inbounds for i in eachindex(x)
        y[i] = clenshaw(a, (x[i] - mid) * invh)
    end
    return y
end

function series_settled(ac, ae, esq)
    scale = 1.0
    @inbounds for v in esq
        av = abs(v)
        av > scale && (scale = av)
    end
    return abs(ac[end]) <= 1e-9 && abs(ae[end]) <= 1e-9 * scale
end

function pair_settled(a::Float64, b::Float64)
    fa = isfinite(a)
    fb = isfinite(b)
    if fa && fb
        m = max(abs(a), abs(b))
        return abs(a - b) <= REL_QUAD * m
    end
    if !fa && !fb
        return (isnan(a) && isnan(b)) || (isinf(a) && isinf(b) && signbit(a) == signbit(b))
    end
    return false
end

function outputs_settled(a, b)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a, b)
        pair_settled(a[i], b[i]) || return false
    end
    return true
end

# Log characteristic function of the centered Gaussian quadratic form, and the
# partial-moment factor β, at a real frequency t. s = i t.
function gauss_phase(L::Layout, t::Float64)
    s = IM * t
    s2 = s * s
    t1 = zero(ComplexF64)
    t2 = zero(ComplexF64)
    t3 = zero(ComplexF64)
    t4 = zero(ComplexF64)
    a2p = zero(ComplexF64)
    a1p = zero(ComplexF64)
    lr = zero(ComplexF64)
    @inbounds for j in eachindex(L.omega)
        nu = 1.0 / (1.0 - 2.0 * L.omega[j] * s)
        nu2 = nu * nu
        t1 += L.d2[j] * nu
        t2 += L.e2[j] * nu
        t3 += L.de[j] * nu
        t4 += log(nu)
        a2p += s * L.d2[j] * nu + s2 * L.d2[j] * L.omega[j] * nu2
        a1p += s * L.e2[j] * nu + s2 * L.e2[j] * L.omega[j] * nu2
        lr += 2.0 * s * L.de[j] * nu + 2.0 * s2 * L.de[j] * L.omega[j] * nu2 + L.omega[j] * nu
    end
    t1 += L.D2z
    t2 += L.E2z
    t3 += L.DEz
    a2p += s * L.D2z
    a1p += s * L.E2z + L.k
    lr += 2.0 * s * L.DEz + L.ccoef
    lm0 = s * L.ccoef + s2 * t3 + 0.5 * t4 + L.k * s + 0.5 * s2 * (t1 + t2)
    beta = a2p + lr
    L.need_a1 && (beta += a1p)
    return lm0, beta
end

function gauss_tail_pair(L::Layout, Ic::Float64, Ip::Float64)
    cval = 0.5 + rp * Ic
    es = (0.5 * L.M20 + rp * Ip) / cval + L.kk
    return cval, es
end

# ∫_T^∞ h(t) exp(-i q t) dt, three terms, h and the partial-moment numerator.
function gauss_ibp!(L::Layout, q::Float64, T::Float64)
    δ = min(1e-7 * (1.0 + T), 0.05 * T)
    function hv(t)
        lm0, beta = gauss_phase(L, t)
        re = real(lm0)
        if re < -700.0 || re > 700.0
            z = zero(ComplexF64)
            return z, z
        end
        e = exp(lm0)
        return e / t, e * beta / t
    end
    hm, pm = hv(T - δ)
    h0, p0 = hv(T)
    hp, pp = hv(T + δ)
    h1 = (hp - hm) / (2.0 * δ)
    p1 = (pp - pm) / (2.0 * δ)
    h2 = (hp - 2.0 * h0 + hm) / (δ * δ)
    p2 = (pp - 2.0 * p0 + pm) / (δ * δ)
    iq = IM * q
    iq2 = iq * iq
    osc = exp(-IM * q * T)
    tc = osc * (h0 / iq + h1 / iq2 + h2 / (iq2 * iq))
    tp = osc * (p0 / iq + p1 / iq2 + p2 / (iq2 * iq))
    return imag(tc), imag(tp)
end

function gauss_panel_osc(L::Layout, q::Float64)
    aq = abs(q)
    T = GAUSS_XT / aq
    wosc = π / (2.0 * aq)
    a = 0.0
    Ic = 0.0
    Ip = 0.0
    @inbounds while a < T
        w = min(wosc, max(1.0, 0.25 * a), T - a)
        w <= 0.0 && break
        mid = a + 0.5 * w
        half = 0.5 * w
        for i in eachindex(GAUSS_GX)
            t = mid + half * GAUSS_GX[i]
            wt = half * GAUSS_GW[i]
            lm0, beta = gauss_phase(L, t)
            lm = lm0 - IM * (t * q)
            Ic += wt * imag_exp(lm) / t
            Ip += wt * imag_exp_mul(lm, beta) / t
        end
        a += w
    end
    tc, tp = gauss_ibp!(L, q, T)
    return gauss_tail_pair(L, Ic + tc, Ip + tp)
end

# q = 0 has no Gil-Pelaez oscillation. The log substitution integrates the
# algebraic tail of a characteristic function that decays as a power.
function gauss_panel_zero(L::Layout)
    Ic = 0.0
    Ip = 0.0
    a = 0.0
    @inbounds while a < 1.0
        w = min(0.05, 1.0 - a)
        w <= 0.0 && break
        mid = a + 0.5 * w
        half = 0.5 * w
        for i in eachindex(GAUSS_GX)
            t = mid + half * GAUSS_GX[i]
            wt = half * GAUSS_GW[i]
            lm0, beta = gauss_phase(L, t)
            Ic += wt * imag_exp(lm0) / t
            Ip += wt * imag_exp_mul(lm0, beta) / t
        end
        a += w
    end
    a = 0.0
    @inbounds while a < 80.0
        w = min(0.5, 80.0 - a)
        w <= 0.0 && break
        mid = a + 0.5 * w
        half = 0.5 * w
        for i in eachindex(GAUSS_GX)
            z = mid + half * GAUSS_GX[i]
            wt = half * GAUSS_GW[i]
            t = exp(z)
            lm0, beta = gauss_phase(L, t)
            Ic += wt * imag_exp(lm0)
            Ip += wt * imag_exp_mul(lm0, beta)
        end
        a += w
    end
    return gauss_tail_pair(L, Ic, Ip)
end

function gauss_kernel_var(L::Layout)
    # sum (d + e)^2 on the kernel. Positive means a Gaussian factor exp(-c t^2).
    return L.D2z + L.E2z + 2.0 * L.DEz
end

function gauss_is_constant(L::Layout)
    return L.kind == KIND_GAUSS && isempty(L.omega) && gauss_kernel_var(L) <= 0.0
end

function gauss_needs_panel(L::Layout)
    L.kind == KIND_GAUSS || return false
    gauss_is_constant(L) && return false
    return gauss_kernel_var(L) <= 0.0 && L.ub > GAUSS_FAST_UB
end

function eval_gauss_constant!(L::Layout, x, ccdf, es)
    @inbounds for qi in eachindex(x)
        q = x[qi] - L.kk
        if q < 0.0
            ccdf[qi] = 1.0
            es[qi] = L.kk
        else
            ccdf[qi] = 0.0
            es[qi] = NaN
        end
    end
    return nothing
end

function eval_gauss_panel!(L::Layout, x, ccdf, es)
    @inbounds for qi in eachindex(x)
        q = x[qi] - L.kk
        if abs(q) < 1e-12
            ccdf[qi], es[qi] = gauss_panel_zero(L)
        else
            ccdf[qi], es[qi] = gauss_panel_osc(L, q)
        end
    end
    return nothing
end

# Compare n nodes with the next refinement on these thresholds.
# Remember the coarse order that agreed, and return the finer values.
# At the cap, return the finest table even if the test is still open.
function certify!(L::Layout, x, ccdf, es)
    if gauss_needs_panel(L)
        eval_gauss_panel!(L, x, ccdf, es)
        return nothing
    end
    n = length(x)
    n <= 0 && return nothing
    nn = L.nlo > 0 ? L.nlo : default_nnode(L.kind)
    nn > NNODE_CAP && (nn = NNODE_CAP)
    nn < 1 && (nn = 1)
    cc_c = Vector{Float64}(undef, n)
    es_c = Vector{Float64}(undef, n)
    cc_f = Vector{Float64}(undef, n)
    es_f = Vector{Float64}(undef, n)
    eval_direct!(materialize(L, nn), x, cc_c, es_c)
    while nn < NNODE_CAP
        n2 = nn * 2
        n2 > NNODE_CAP && (n2 = NNODE_CAP)
        eval_direct!(materialize(L, n2), x, cc_f, es_f)
        if outputs_settled(cc_c, cc_f) && outputs_settled(es_c, es_f)
            L.nlo = nn
            copyto!(ccdf, cc_f)
            copyto!(es, es_f)
            return nothing
        end
        if n2 == NNODE_CAP
            copyto!(ccdf, cc_f)
            copyto!(es, es_f)
            return nothing
        end
        cc_c, cc_f = cc_f, cc_c
        es_c, es_f = es_f, es_c
        nn = n2
    end
    copyto!(ccdf, cc_c)
    copyto!(es, es_c)
    return nothing
end

function eval_quad!(L::Layout, x, ccdf, es)
    n = length(x)
    n <= 0 && return nothing
    if gauss_is_constant(L)
        eval_gauss_constant!(L, x, ccdf, es)
        return nothing
    end
    if n <= 24
        certify!(L, x, ccdf, es)
        return nothing
    end
    xmin = xmax = x[1]
    @inbounds for i = 2:n
        x[i] < xmin && (xmin = x[i])
        x[i] > xmax && (xmax = x[i])
    end
    if xmax - xmin <= 1e-14 * (1.0 + abs(xmax))
        one = Vector{Float64}(undef, 1)
        one[1] = x[1]
        cc = Vector{Float64}(undef, 1)
        ee = Vector{Float64}(undef, 1)
        certify!(L, one, cc, ee)
        fill!(ccdf, cc[1])
        fill!(es, ee[1])
        return nothing
    end
    mid = 0.5 * (xmin + xmax)
    half = 0.5 * (xmax - xmin)
    # A smooth tail is analytic in the threshold, so a short Chebyshev grid
    # reproduces the node rule. A sharp bend is not, and is integrated directly.
    # The node count is certified on the abscissae of this call, and the coarse
    # order is kept for the next degree so a hard tail is not restarted at 32.
    m = min(cheb_order(L.kind), n)
    while true
        xq = Vector{Float64}(undef, m)
        @inbounds for j = 1:m
            xq[j] = mid + half * cos(π * (j - 1 + 0.5) / m)
        end
        ccq = Vector{Float64}(undef, m)
        esq = Vector{Float64}(undef, m)
        certify!(L, xq, ccq, esq)
        ac = Vector{Float64}(undef, m)
        ae = Vector{Float64}(undef, m)
        cheb_coeffs!(ac, ccq)
        cheb_coeffs!(ae, esq)
        if series_settled(ac, ae, esq)
            interp!(ccdf, x, ccq, mid, half)
            interp!(es, x, esq, mid, half)
            return nothing
        end
        if m >= n || m >= 96
            certify!(L, x, ccdf, es)
            return nothing
        end
        nxt = min(n, 96, m + max(8, m ÷ 2))
        if nxt <= m
            certify!(L, x, ccdf, es)
            return nothing
        end
        m = nxt
    end
end

function quadrature_eval(x::AbstractVector, a0, a, A, C, mu, gam, lam, chi, psi)
    xv = collect(Float64, x)
    L = prepare(a0, a, A, C, mu, gam, lam, chi, psi)
    ccdf = similar(xv)
    es = similar(xv)
    eval_quad!(L, xv, ccdf, es)
    return ccdf, es
end

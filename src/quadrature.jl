# Mapped Gauss–Legendre evaluation of the Gil-Pelaez integral.
# The node counts, the truncation, and the Chebyshev grid match es4mgh.c.

const KIND_GENERAL = 0
const KIND_NIG = 1
const KIND_HALF = 2
const KIND_PSI0 = 3
const IM = ComplexF64(0.0, 1.0)
const LOG2F = log(2.0)

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
    (kind == KIND_NIG || kind == KIND_GENERAL) && return 32
    kind == KIND_HALF && return 48
    return 64
end

function cheb_order(kind::Int)
    # First degree tried. The normal-inverse-Gaussian tail has settled by 20.
    # Other laws start higher and still grow until the last coefficient is small.
    kind == KIND_NIG && return 20
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
    nn = default_nnode(kind)
    ub = integration_ub(omega_all)
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
            nu = 1.0 / (1.0 - 2.0 * omega[j] * s)
            nu2 = nu * nu
            t1 += d2[j] * nu
            t2 += e2[j] * nu
            t3 += de[j] * nu
            t4 += log(nu)
            a2p += s * d2[j] * nu + s2 * d2[j] * omega[j] * nu2
            a1p += s * e2[j] * nu + s2 * e2[j] * omega[j] * nu2
            lr += 2.0 * s * de[j] * nu + 2.0 * s2 * de[j] * omega[j] * nu2 + omega[j] * nu
        end
        t1 += D2z
        t2 += E2z
        t3 += DEz
        a2p += s * D2z
        a1p += s * E2z + k
        lr += 2.0 * s * DEz + ccoef
        chi_base[i] = chi - s2 * t1
        pnode = psi - 2.0 * (k * s + 0.5 * s2 * t2)
        psi_node[i] = pnode
        lrho[i] = s * ccoef + s2 * t3 + 0.5 * t4
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
    lm1 = lklam_real_f(lam + 1.0, chi, psi) - LK2
    lm2 = lklam_real_f(lam + 2.0, chi, psi) - LK2
    sum_om = isempty(omega_all) ? 0.0 : sum(omega_all)
    # The skewness term carries E[W^2]. When that coefficient is zero the
    # product is zero even if the moment is infinite (integer λ, ψ = 0).
    skew = k == 0.0 ? 0.0 : real(exp(lm2) * k)
    M20 = skew + real(exp(lm1) * (ccoef + sum_om))
    return Prepared(kind, kk, lam, LK2, M20, need_a1, nn, u, wo, chi_base, psi_node, lrho, a2pv, a1pv, lrp, log_psi)
end

function finish_tail(E::Prepared, Ic::Float64, Ip::Float64)
    cval = 0.5 + rp * Ic
    es = (0.5 * E.M20 + rp * Ip) / cval + E.kk
    return cval, es
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
    if E.kind == KIND_NIG
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
    if E.kind == KIND_NIG
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

function eval_quad!(E::Prepared, x, ccdf, es)
    n = length(x)
    n <= 0 && return nothing
    if n <= 24
        eval_direct!(E, x, ccdf, es)
        return nothing
    end
    xmin = xmax = x[1]
    @inbounds for i = 2:n
        x[i] < xmin && (xmin = x[i])
        x[i] > xmax && (xmax = x[i])
    end
    if xmax - xmin <= 1e-14 * (1.0 + abs(xmax))
        cc, ee = eval_point(E, x[1] - E.kk)
        fill!(ccdf, cc)
        fill!(es, ee)
        return nothing
    end
    mid = 0.5 * (xmin + xmax)
    half = 0.5 * (xmax - xmin)
    # A smooth tail is analytic in the threshold, so a short Chebyshev grid
    # reproduces the node rule. A sharp bend is not, and is integrated directly.
    m = min(cheb_order(E.kind), n)
    while true
        xq = Vector{Float64}(undef, m)
        @inbounds for j = 1:m
            xq[j] = mid + half * cos(π * (j - 1 + 0.5) / m)
        end
        ccq = Vector{Float64}(undef, m)
        esq = Vector{Float64}(undef, m)
        eval_direct!(E, xq, ccq, esq)
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
            eval_direct!(E, x, ccdf, es)
            return nothing
        end
        nxt = min(n, 96, m + max(8, m ÷ 2))
        if nxt <= m
            eval_direct!(E, x, ccdf, es)
            return nothing
        end
        m = nxt
    end
end

function quadrature_eval(x::AbstractVector, a0, a, A, C, mu, gam, lam, chi, psi)
    xv = collect(Float64, x)
    E = prepare(a0, a, A, C, mu, gam, lam, chi, psi)
    ccdf = similar(xv)
    es = similar(xv)
    eval_quad!(E, xv, ccdf, es)
    return ccdf, es
end

# Values pinned here are from QuadraticFormsMGHyp 0.1.1 (the C kernel),
# except the psi = 0 shortfall, which that kernel returns as NaN.
using LinearAlgebra
using QuadraticFormsMGHyp

const ATOL = 1e-9
const ES_ATOL = 1e-8

function twosls_case(nu, b)
    n = 25
    strength = 0.5
    instrument = zeros(n)
    instrument[1] = sqrt(strength)
    projection = zeros(n, n)
    projection[1, 1] = 1.0
    eye = Matrix{Float64}(I, n, n)
    cov = ((nu - 2.0) / nu) * ones(2, 2)
    symmetric = 0.5 * (cov + transpose(cov))
    F = eigen(Symmetric(symmetric))
    scaled = F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
    factor = scaled * transpose(F.vectors)
    scale = kron(factor, eye)
    a_vec = vcat(instrument, -2.0 * b * instrument)
    quad = [zeros(n, n) 0.5 * projection; 0.5 * projection (-b * projection)]
    return qfmgh(0.0, -b * strength, a_vec, quad, scale, zeros(2n), zeros(2n), -0.5 * nu, nu, 0.0)
end

@testset "mapped quadrature" begin
    z = zeros(10)
    A = 0.5 * Matrix{Float64}(I, 10, 10)
    C = Matrix{Float64}(I, 10, 10)
    x = range(3.5, 17.5; length = 1401)
    c1, e1 = qfmgh(3.5, 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    cmid, emid = qfmgh(x[701], 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    cend, eend = qfmgh(17.5, 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    @test c1 isa Float64
    @test c1 ≈ 0.45063991441010837 atol = ATOL
    @test e1 ≈ 8.90139728539701 atol = ES_ATOL
    ccdf, es = qfmgh(x, 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    @test ccdf[1] ≈ c1 atol = 1e-8
    @test ccdf[701] ≈ cmid atol = 1e-8
    @test ccdf[end] ≈ cend atol = 1e-8
    @test es[1] ≈ e1 atol = 1e-7
    @test es[701] ≈ emid atol = 1e-7
    @test es[end] ≈ eend atol = 1e-7
    @test issorted(reverse(ccdf))

    # Twenty-four thresholds are integrated directly. Twenty-five use the Chebyshev grid.
    c24, e24 = qfmgh(range(3.5, 17.5; length = 24), 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    @test c24[1] ≈ 0.45063991441010837 atol = ATOL
    @test c24[end] ≈ 0.03923559352643207 atol = ATOL
    @test e24[1] ≈ 8.90139728539701 atol = ES_ATOL
    c25, e25 = qfmgh(range(3.5, 17.5; length = 25), 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    @test c25[1] ≈ ccdf[1] atol = 1e-8
    @test e25[1] ≈ es[1] atol = 1e-7
    @test c25[end] ≈ ccdf[end] atol = 1e-8

    flat, flates = qfmgh(fill(6.0, 40), 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    one, onees = qfmgh(6.0, 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
    @test all(v -> v ≈ one, flat)
    @test all(v -> v ≈ onees, flates)
    @test one ≈ 0.25373837265745125 atol = ATOL
    @test onees ≈ 12.24495342880734 atol = ES_ATOL

    # Chebyshev against the same rule evaluated at the threshold itself.
    xs = collect(x)
    for i in (1, 50, 400, 700, 1000, 1401)
        ci, ei = qfmgh(xs[i], 0.0, z, A, C, z, z, -0.5, 1.0, 1.0)
        @test ccdf[i] ≈ ci atol = 1e-7
        @test es[i] ≈ ei atol = 1e-6
    end

    gam = collect(range(-0.2, 0.3; length = 10))
    mu = collect(range(0.0, 0.1; length = 10))
    a = collect(range(-0.1, 0.2; length = 10))
    xs = range(0.0, 12.0; length = 80)
    cs1, es1 = qfmgh(0.0, 0.4, a, A, C, mu, gam, -0.5, 1.2, 0.8)
    csend, esend = qfmgh(12.0, 0.4, a, A, C, mu, gam, -0.5, 1.2, 0.8)
    @test cs1 ≈ 0.9999999966436798 atol = ATOL
    @test es1 ≈ 7.331363902541845 atol = ES_ATOL
    cs, ess = qfmgh(xs, 0.4, a, A, C, mu, gam, -0.5, 1.2, 0.8)
    @test cs[1] ≈ cs1 atol = 1e-8
    @test cs[end] ≈ csend atol = 1e-8
    @test ess[1] ≈ es1 atol = 1e-7
    @test ess[end] ≈ esend atol = 1e-7

    # Scalar calls use the node rule itself. The digits are that rule.
    # A long grid must stay within 1e-8 of it. A fixed six-node interpolant does not.
    ch1, eh1 = qfmgh(1.0, 0.0, z, A, C, z, z, -1.5, 2.0, 1.5)
    ch20, eh20 = qfmgh(20.0, 0.0, z, A, C, z, z, -1.5, 2.0, 1.5)
    @test ch1 ≈ 0.8939020478983136 atol = ATOL
    @test eh1 ≈ 4.00912199544807 atol = ES_ATOL
    ch, eh = qfmgh(range(1.0, 20.0; length = 200), 0.0, z, A, C, z, z, -1.5, 2.0, 1.5)
    @test ch[1] ≈ ch1 atol = 1e-8
    @test ch[end] ≈ ch20 atol = 1e-8
    @test eh[1] ≈ eh1 atol = 1e-7
    @test eh[end] ≈ eh20 atol = 1e-7

    cg1, eg1 = qfmgh(1.0, 0.0, z, A, C, z, z, -1.3, 2.0, 1.1)
    cg20, eg20 = qfmgh(20.0, 0.0, z, A, C, z, z, -1.3, 2.0, 1.1)
    @test cg1 ≈ 0.9186095325658532 atol = ATOL
    @test eg1 ≈ 4.776914154468738 atol = ES_ATOL
    cg, eg = qfmgh(range(1.0, 20.0; length = 200), 0.0, z, A, C, z, z, -1.3, 2.0, 1.1)
    @test cg[1] ≈ cg1 atol = 1e-8
    @test cg[end] ≈ cg20 atol = 1e-8
    @test eg[1] ≈ eg1 atol = 1e-7
    @test eg[end] ≈ eg20 atol = 1e-7

    ci1, ei1 = qfmgh(1.0, 0.0, z, A, C, z, z, -2.0, 3.0, 1.0)
    ci15, ei15 = qfmgh(15.0, 0.0, z, A, C, z, z, -2.0, 3.0, 1.0)
    @test ci1 ≈ 0.9417926650711783 atol = ATOL
    @test ei1 ≈ 4.674253359582102 atol = ES_ATOL
    ci, ei = qfmgh(range(1.0, 15.0; length = 40), 0.0, z, A, C, z, z, -2.0, 3.0, 1.0)
    @test ci[1] ≈ ci1 atol = 1e-8
    @test ci[end] ≈ ci15 atol = 1e-8
    @test ei[1] ≈ ei1 atol = 1e-7
    @test ei[end] ≈ ei15 atol = 1e-7

    # λ = -2 and ψ = 0: E[W^2] is infinite, but the skewness coefficient is zero.
    # The tail bends sharply, so a long grid is integrated point by point.
    A4 = diagm([0.4, 0.0, 0.2, 0.7])
    p4 = (0.1, [0.2, -0.1, 0.0, 0.3], A4, Matrix{Float64}(I, 4, 4), zeros(4), zeros(4), -2.0, 4.0, 0.0)
    cp1, ep1 = qfmgh(-1.0, p4...)
    cp8, ep8 = qfmgh(8.0, p4...)
    @test cp1 ≈ 0.9999998269559728 atol = ATOL
    @test ep1 ≈ 2.700000681379396 atol = ES_ATOL
    cp, ep = qfmgh(range(-1.0, 8.0; length = 50), p4...)
    @test cp[1] ≈ cp1 atol = ATOL
    @test cp[end] ≈ cp8 atol = ATOL
    @test all(isfinite, ep)
    @test ep[1] ≈ ep1 atol = ES_ATOL
    @test ep[end] ≈ ep8 atol = ES_ATOL
    @test ep1 > 0

    # Portfolio greeks from the Python examples, diagonal and therefore stable.
    a0 = -0.3284225867331385
    avec = [-0.41141088640242746, -0.41141088640242746, -0.41141088640242746,
        -0.41141088640242746, -0.41141088640242746, 0.5885891135975725,
        0.5885891135975725, 0.5885891135975725, 0.5885891135975725, 0.5885891135975725]
    Ad = diagm(fill(0.0091703580324228, 10))
    Cd = diagm(fill(1.8898223650461363, 10))
    cp1, ep1p = qfmgh(3.5, a0, avec, Ad, Cd, z, z, -0.5, 1.0, 1.0)
    @test cp1 ≈ 0.10147755602526368 atol = 1e-8
    @test ep1p ≈ 5.801375708121332 atol = 1e-7
    cport, eport = qfmgh(range(3.5, 17.5; length = 1401), a0, avec, Ad, Cd, z, z, -0.5, 1.0, 1.0)
    @test cport[1] ≈ cp1 atol = 1e-8
    @test eport[1] ≈ ep1p atol = 1e-7
    @test issorted(reverse(cport))

    # Fifty-dimensional 2SLS. Eigenvalues move between BLAS libraries, so the
    # tolerance is the gap observed against Accelerate, with room to spare.
    c0, e0 = twosls_case(3.0, 0.0)
    @test c0 ≈ 0.6540340046251787 atol = 1e-8
    @test e0 ≈ 1.5735661995338892 atol = 1e-8
    c3, e3 = twosls_case(3.0, 3.0)
    @test c3 ≈ 0.07228081754479557 atol = 1e-5
    @test e3 ≈ 0.041358477407942074 atol = 1e-5
    c9, e9 = twosls_case(9.0, 0.0)
    @test c9 ≈ 0.7216659250848524 atol = 1e-8
    @test e9 ≈ 1.4181174952966877 atol = 1e-8
end

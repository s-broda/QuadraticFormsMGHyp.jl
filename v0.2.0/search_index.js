var documenterSearchIndex = {"docs":
[{"location":"#QuadraticFormsMGHyp.jl-1","page":"Home","title":"QuadraticFormsMGHyp.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [QuadraticFormsMGHyp]","category":"page"},{"location":"#QuadraticFormsMGHyp.qfmgh","page":"Home","title":"QuadraticFormsMGHyp.qfmgh","text":"qfmgh(x, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)\n\nSurvivor function P(L>x) = 1-F(x) and tail conditional mean, E(L|L>x), of\n\nL = a0 + a' * X + X' * A * X, where:\n\nX = mu + W * gam + sqrt(W) * C * Z,    Z ~ N(0, I), and W ~ GIG(lam, chi, psi), i.e.,    X is distributed as multivariate GHyp.\n\nKeyword arguments:\n\n`do_spa`: whether to return the exact result or a saddlepoint approximation\n`order`: order of the saddlepoint approximation\n\n(c) 2020 S.A. Broda\n\n\n\n\n\n","category":"function"}]
}
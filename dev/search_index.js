var documenterSearchIndex = {"docs":
[{"location":"#QuadraticFormsMGHyp.jl-1","page":"Home","title":"QuadraticFormsMGHyp.jl","text":"","category":"section"},{"location":"#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package implements the algorithms from our paper On Quadratic Forms in Multivariate Generalized Hyperbolic Random Vectors, which deals with tail probabilities and partial moments of quadratic forms.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Consider the random variable","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Lequiv a_0+mathbfa^mathrmscriptscriptstyle TX+X^mathrmscriptscriptstyle TmathbfAX","category":"page"},{"location":"#","page":"Home","title":"Home","text":"a quadratic plus a linear form in the random vector X. Suppose Xsim mathrmMGHyp(boldsymbolmumathbfCboldsymbolgammalambdachipsi); that is, X has a d-variate generalized hyperbolic distribution with stochastic representation","category":"page"},{"location":"#","page":"Home","title":"Home","text":"X=boldsymbolmu+Y boldsymbolgamma +surdYmathbfCZ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where Z has a d-variate standard Normal distribution, boldsymbolmu and  boldsymbolgamma are constant d-vectors, mathbfC is a dtimes d matrix, and Y has a univariate generalized inverse Gaussian distribution with density","category":"page"},{"location":"#","page":"Home","title":"Home","text":"f_GIG(ylambdachipsi)equivfracy^lambda-1k_lambda(chipsi)expleft-frac12left(chi y^-1+psi yright)right","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where","category":"page"},{"location":"#","page":"Home","title":"Home","text":"k_lambda(chipsi)equivbegincasesfracpsi2^-lambdaGamma(lambda)text if chi=0\nfracchi2^lambdaGamma(-lambda)text if psi=0\n2left(fracchipsiright)^lambda2K_lambda(sqrtchipsi) text if chineq0 text and psineq0endcases","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Here, K_lambda(z) is the modified Bessel function of the second kind of order nu.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The package exports a single function, qfmgh. Its signature is","category":"page"},{"location":"#","page":"Home","title":"Home","text":"qfmgh(x, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The keyword argument do_spa controls whether an exact result or a saddlepoint approximation is computed. The order of the latter is controlled with the second keyword argument, order, which can be either 1 or 2.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The function returns a touple containing the tail probability mathbbPLx and the tail conditional mean mathbbELmid Lx .","category":"page"},{"location":"#Docstrings-1","page":"Home","title":"Docstrings","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Modules = [QuadraticFormsMGHyp]","category":"page"},{"location":"#QuadraticFormsMGHyp.qfmgh","page":"Home","title":"QuadraticFormsMGHyp.qfmgh","text":"qfmgh(x, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)\n\nSurvivor function P(L>x) = 1-F(x) and tail conditional mean, E(L|L>x), of\n\nL = a0 + a' * X + X' * A * X, where:\n\nX = mu + W * gam + sqrt(W) * C * Z,    Z ~ N(0, I), and W ~ GIG(lam, chi, psi), i.e.,    X is distributed as multivariate GHyp.\n\nKeyword arguments:\n\n`do_spa`: whether to return the exact result or a saddlepoint approximation\n`order`: order of the saddlepoint approximation\n\n(c) 2020 S.A. Broda\n\n\n\n\n\n","category":"function"}]
}

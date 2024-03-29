# QuadraticFormsMGHyp.jl
## Introduction
This package implements the algorithms from our paper [On Quadratic Forms in Multivariate Generalized Hyperbolic Random Vectors](https://dx.doi.org/10.2139/ssrn.3369208), which deals with tail probabilities and partial moments of quadratic forms. The algorithms generalize those of Imhof (Biometrika, 1961) and Broda (Mathematical Finance, 2012).

Consider the random variable
```math
L\equiv a_0+\mathbf{a}^{\mathrm{\scriptscriptstyle T}}X+X^{\mathrm{\scriptscriptstyle T}}\mathbf{A}X,
```
a quadratic plus a linear form in the random vector ``X\sim \mathrm{MGHyp}(\boldsymbol{\mu},\mathbf{C},\boldsymbol{\gamma},\lambda,\chi,\psi)``; i.e., ``X`` has a ``d``-variate generalized hyperbolic distribution with stochastic representation
```math
X=\boldsymbol{\mu}+Y \boldsymbol{\gamma} +\surd{Y}\mathbf{C}Z,
```
where ``Z`` has a ``d``-variate standard Normal distribution, ``\boldsymbol{\mu}`` and  ``\boldsymbol{\gamma}`` are constant ``d``-vectors, ``\mathbf{C}`` is a ``d\times d`` matrix, and ``Y`` has a univariate generalized inverse Gaussian distribution with density
```math
f_{GIG}(y;\lambda,\chi,\psi)\equiv\frac{y^{\lambda-1}}{k_\lambda(\chi,\psi)}\exp\left\{-\frac{1}{2}\left(\chi y^{-1}+\psi y\right)\right\},
```
where
```math
k_\lambda(\chi,\psi)\equiv\begin{cases}\frac{\psi}{2}^{-\lambda}\Gamma(\lambda),\text{ if }\chi=0\\
\frac{\chi}{2}^{\lambda}\Gamma(-\lambda),\text{ if }\psi=0\\
2\left(\frac{\chi}{\psi}\right)^{\lambda/2}K_\lambda(\sqrt{\chi\psi}), \text{ if }\chi\neq0 \text{ and }\psi\neq0,\end{cases}
```
and ``K_\lambda(z)`` is the modified Bessel function of the second kind of order ``\nu``.

The generalized hyperbolic distribution contains as special cases, among others, the Variance-Gamma (``\lambda>0``), Student's ``t`` (``\lambda=-\nu/2, \chi=\nu, \psi=0``), Normal Inverse Gaussian (``\lambda=-1/2``), and Hyperbolic (``\lambda=1``) distributions.

## Installation
`QuadraticFormsMGHyp` is a registered Julia package. It can be installed with `using Pkg; Pkg.add("QuadraticFormsMGHyp")`. Figure 1 in the paper, reproduced in the logo above, can be generated with `using Pkg; Pkg.test("QuadraticFormsMGHyp")`.

## Usage
The package exports a single function, [`qfmgh`](@ref). Its signature is

```julia
qfmgh(x, a0, a, A, C, mu, gam, lam, chi, psi; do_spa=false, order=2)
```
The first argument, `x`, can be passed as either a scalar or a vector. The latter is more efficient than the more Julian way of calling `qfmgh.(xvec, ...)`, because certain parts of the computation, such as calculating the eigenvalues of ``A``, can be hoisted out of the loop over `xvec`.

The keyword argument `do_spa` controls whether an exact result or a saddlepoint approximation is computed. The order of the latter is controlled with the second keyword argument, `order`, which can be either 1 or 2.

The function returns a tuple containing (vectors of) the tail probability ``\mathbb{P}[L>x]``, and the tail conditional mean ``\mathbb{E}[L\mid L>x ]``.

## Docstrings
```@autodocs
Modules = [QuadraticFormsMGHyp]
```

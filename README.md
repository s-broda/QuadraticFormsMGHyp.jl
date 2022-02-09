# QuadraticFormsMGHyp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://s-broda.github.io/QuadraticFormsMGHyp.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://s-broda.github.io/QuadraticFormsMGHyp.jl/dev)
[![Build Status](https://travis-ci.org/s-broda/QuadraticFormsMGHyp.jl.svg?branch=master)](https://travis-ci.org/s-broda/QuadraticFormsMGHyp.jl)
[![Build status (Linux, MacOS)](https://github.com/s-broda/ARCHModels.jl/workflows/CI/badge.svg)](https://github.com/s-broda/ARCHModels.jl/actions?query=workflow%3ACI)
[![Build Status (Windows)](https://ci.appveyor.com/api/projects/status/github/s-broda/QuadraticFormsMGHyp.jl?svg=true)](https://ci.appveyor.com/project/s-broda/QuadraticFormsMGHyp-jl)
[![Codecov](https://codecov.io/gh/s-broda/QuadraticFormsMGHyp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/s-broda/QuadraticFormsMGHyp.jl)
[![pkgeval](https://juliahub.com/docs/QuadraticFormsMGHyp/pkgeval.svg)](https://juliahub.com/ui/Packages/QuadraticFormsMGHyp/vxdaX)

A package for evaluating tail probabilities and partial moments for random vectors in multivariate generalized hyperbolic random vectors. Matlab and Fortran code is available [here](https://github.com/s-broda/es4mgh).

# Installation
The package can be installed with `using Pkg; Pkg.add("QuadraticFormsMGHyp")`.

# Citation
If you use this package in your research, then please consider citing [our paper](https://doi.org/10.1093/biomet/asaa067). The figures in the paper can be recreated by running `using Pkg; Pkg.test("QuadraticFormsMGHyp")`.

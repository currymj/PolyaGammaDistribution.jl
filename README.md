# PolyaGammaDistribution

[![Build Status](https://travis-ci.org/currymj/PolyaGammaDistribution.jl.svg?branch=master)](https://travis-ci.org/currymj/PolyaGammaDistribution.jl)

[![Coverage Status](https://coveralls.io/repos/currymj/PolyaGammaDistribution.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/currymj/PolyaGammaDistribution.jl?branch=master)

[![codecov.io](http://codecov.io/github/currymj/PolyaGammaDistribution.jl/coverage.svg?branch=master)](http://codecov.io/github/currymj/PolyaGammaDistribution.jl?branch=master)

[![Build status](https://ci.appveyor.com/api/projects/status/1cbrjw3gotlvps31?svg=true)](https://ci.appveyor.com/project/currymj/polyagammadistribution-jl)

# What is this?

This repository is still under active development, but when it's done it will provide tools for sampling from and computing moments of the Pólya-Gamma distribution, as described in [this paper by Polson et al](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001). They provide the [BayesLogit R package](https://cran.r-project.org/web/packages/BayesLogit/index.html) which my code at the moment essentially copies from, hence the GPL license.

It uses the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl.git) interface.

As far as I can tell, it works correctly right now, for distributions PG(b, c) with integer-valued b parameter. Work remains to ensure this is true, and to optimize the code and make it more Julian and less of a literal
translation of the R code.

The integer-valued case is what's needed for almost all of the statistical applications. It is possible to generalize, less efficiently, and this may be implemented as well.

# How to use it

```julia
using PolyaGammaDistribution

pgdist = PolyaGamma(1,3.0)

singlesample = rand(pgdist)
manysamples = rand(pgdist, 10)
```

# What is the point of this?

Suppose we have some function, often some sort of distance between two points, that outputs a real value. We'd like a probability, so we push it through a function that squashes from 0 to 1. Often the logistic aka sigmoid function is used, as in logistic regression or the common activations on the end of a neural network. This is an easy thing to do, although it usually ends up requiring some kind of numerical method to maximize the output.

However, if you have a model where you'd like to use a sampler to approximate a posterior distribution (rather than just maximizing) then it becomes much trickier. The Polson et al paper lays out an efficient scheme for sampling in these situations by adding some Pólya-Gamma distributed latent variables.

# Acknowledgments, etc.

Apologies to George Pólya for misspelling his name, but it seemed like putting non-ASCII characters in a package/module name was asking for trouble down the line.

N. G. Polson, J. G. Scott, and J. Windle, “Bayesian Inference for Logistic Models Using Pólya–Gamma Latent Variables,” Journal of the American Statistical Association, vol. 108, no. 504, pp. 1339–1349, Dec. 2013.

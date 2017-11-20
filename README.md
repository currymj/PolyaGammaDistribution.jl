# PolyaGammaDistribution

[![Build Status](https://travis-ci.org/currymj/PolyaGamma.jl.svg?branch=master)](https://travis-ci.org/currymj/PolyaGamma.jl)

[![Coverage Status](https://coveralls.io/repos/currymj/PolyaGamma.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/currymj/PolyaGamma.jl?branch=master)

[![codecov.io](http://codecov.io/github/currymj/PolyaGamma.jl/coverage.svg?branch=master)](http://codecov.io/github/currymj/PolyaGamma.jl?branch=master)

# What is this?

This repository is still under rapid and active development, but when it's done it will provide tools for sampling from and computing moments of the Pólya-Gamma distribution, as described in [this paper by Polson et al](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001). They provide the [BayesLogit R package](https://cran.r-project.org/web/packages/BayesLogit/index.html) which my code at the moment essentially copies from, hence the GPL license.


# What is the point of this?

Suppose we have some function, often some sort of distance between two points, that outputs a real value. We'd like a probability, so we push it through a function that squashes from 0 to 1. Often the logistic aka sigmoid function is used, as in logistic regression or the common activations on the end of a neural network. This is an easy thing to do, although it usually ends up requiring some kind of numerical method to maximize the output.

However, if you have a model where you'd like to use a sampler to approximate a posterior distribution (rather than just maximizing) then it becomes much trickier. The Polson et al paper lays out an efficient scheme for sampling in these situations by adding some Pólya-Gamma distributed latent variables.

# Acknowledgments, etc.

Apologies to George Pólya for misspelling his name, but it seemed like putting non-ASCII characters in a package/module name was asking for trouble down the line.

N. G. Polson, J. G. Scott, and J. Windle, “Bayesian Inference for Logistic Models Using Pólya–Gamma Latent Variables,” Journal of the American Statistical Association, vol. 108, no. 504, pp. 1339–1349, Dec. 2013.


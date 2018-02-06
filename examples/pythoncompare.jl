# this is just a simple script i wrote to compare outputs to pypolyagamma with KS test
using PyCall
using HypothesisTests
using PolyaGammaDistribution

@pyimport pypolyagamma

function python_samples(N)
    pysampler = pypolyagamma.PyPolyaGamma()

    [pysampler[:pgdraw](1,3) for _=1:N]
end

function julia_samples(N)
    juliasampler = PolyaGamma(1.0, 3.0)
    rand(juliasampler, N)
end


function kscompare(N)
    p = python_samples(N)
    j = julia_samples(N)
    ApproximateTwoSampleKSTest(j, p)
end

test = kscompare(10000)

pvalue(test)

# the null hypothesis is that they are both from the same distribution,
# so very large p values is good.

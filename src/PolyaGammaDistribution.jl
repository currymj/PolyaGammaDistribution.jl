module PolyaGammaDistribution
using Distributions

const TRUNC = 0.64
const cutoff = 1 / TRUNC

struct PolyaGamma{T<:Real} <: ContinuousUnivariateDistribution
    b::T
    c::T
end

function Distributions.mean(d::PolyaGamma)
    (d.b / (2.0*d.c)) * tanh(d.c / 2.0)
end

# functions from BayesLogit

# cdf of Inverse Gaussian, already helpfully given to us
function pigauss(x, μ, λ)
    cdf(InverseGaussian(μ, λ), x)
end

function rtigauss(z, r=TRUNC)
    tdist = Truncated(InverseGaussian(1/abs(z), 1.0), 0.0, r)
    rand(tdist)
end

function rand_gammasum(d::PolyaGamma{T}, num=1, n=1, z=0.0, trunc=200) where T<:Real
    ci = ((1:num).-(1/2)).^2 * pi^2 * 4
    ai = ci + z.^2
    w = zeros(T, num)
    for i=1:num
        w[i] = 2.0 * sum(rand(Gamma(n),trunc)./ai)
    end
    w
end

export PolyaGamma, mean
end # module

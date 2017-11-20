module PolyaGamma
using Distributions

struct PolyaGammaDistribution{T<:Real} <: ContinuousUnivariateDistribution
    b::T
    c::T
end

function mean(d::PolyaGammaDistribution)
    (d.b / (2.0*d.c)) * tanh(d.c / 2.0)
end
# package code goes here
export PolyaGammaDistribution, mean
end # module

using PolyaGammaDistribution
using Distributions
using Plots

# copy example from pypolyagamma

N = 10
mu = 0.0
sigmasq = 1.0

σ(x) = 1.0/(1 + exp(-x))
xtrue = rand(Normal(mu, sqrt(sigmasq)))
ptrue = σ(xtrue)

y = rand(Binomial(N, ptrue))
N_samples = 10000
xs = zeros(N_samples)
ωs = ones(N_samples)

for i=2:N_samples
    ωs[i] = rand(PolyaGamma(N, xs[i-1]))

    sigmasq_hat = 1./(1. / sigmasq + ωs[i])
    mu_hat = sigmasq_hat * (mu / sigmasq + (y - N / 2.))
    xs[i] = rand(Normal(mu_hat, sqrt(sigmasq_hat)))
end

histogram(xs)

xtrue

mean(xs)
mode(xs)

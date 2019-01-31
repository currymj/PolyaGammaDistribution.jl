using PolyaGammaDistribution
using Test
using Distributions: mean, var
using Random: GLOBAL_RNG

# Distributions.jl is willing to have tests that fail with some small probability
# even when the code is working, so let's just do that here.
@testset begin
    @test mean(PolyaGamma(2, 3.0)) ≈ 0.3017160845482888
    @test mean(PolyaGamma(1, 1.0)) ≈ 0.23105857863000487
    @test var(PolyaGamma(1, 1.0)) ≈ 0.03444664538852302
end

@testset begin
    tol = 1e-3
    @test abs(mean([PolyaGammaDistribution.rtigauss(GLOBAL_RNG, 1.0) for _ in 1:10000]) - .372498) < .005
    @test abs(PolyaGammaDistribution.mass_texpon(0.0) - 0.5776972) < tol
    @test abs(PolyaGammaDistribution.mass_texpon(1.0) - 0.4605903) < tol
    @test abs(PolyaGammaDistribution.mass_texpon(2.0) - 0.2305365) < tol
end

@testset begin
    sigma_tol = 5
    d = PolyaGamma(1, 1.0)
    analytic_mean = mean(d)
    analytic_var = var(d)
    nsamples = 10^6
    # standard error of the sample mean
    standard_err = sqrt(analytic_var / nsamples)
    # standard error of the sample variance
    se_variance = sqrt(2*analytic_var^2/nsamples)
    @test abs(analytic_mean - mean(rand(d, nsamples))) < sigma_tol*standard_err
    @test abs(analytic_var -  var(rand(d, nsamples))) < sigma_tol*se_variance
end

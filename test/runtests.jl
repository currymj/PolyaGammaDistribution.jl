using PolyaGammaDistribution
using Base.Test

# Distributions.jl is willing to have tests that fail with some small probability
# even when the code is working, so let's just do that here.
@testset begin
    @test mean(PolyaGamma(2.0, 3.0)) ≈ 0.3017160845482888
    @test mean(PolyaGamma(1.0, 1.0)) ≈ 0.23105857863000487
    @test var(PolyaGamma(1.0, 1.0)) ≈ 0.03444664538852302
end

@testset begin
    @test abs(mean([PolyaGammaDistribution.rtigauss(1.0) for _ in 1:10000]) - .372498) < .005
end

@testset begin
    const percent_tol = 0.01
    analytic_mean = mean(PolyaGamma(1.0, 1.0))
    analytic_var = var(PolyaGamma(1.0, 1.0))
    @test abs(analytic_mean - mean(rand(PolyaGamma(1.0, 1.0), 1000000))) < abs(percent_tol * analytic_mean)
    @test abs(analytic_var - var(rand(PolyaGamma(1.0, 1.0), 1000000))) < abs(percent_tol * analytic_var)
end
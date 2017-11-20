using PolyaGammaDistribution
using Base.Test

@testset begin
    @test PolyaGammaDistribution.mean(PolyaGammaDistribution.PolyaGamma(2.0, 3.0)) ≈ 0.3017160845482888
    @test PolyaGammaDistribution.mean(PolyaGammaDistribution.PolyaGamma(1.0, 1.0)) ≈ 0.23105857863000487
end
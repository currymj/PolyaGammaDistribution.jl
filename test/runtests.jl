using PolyaGamma
using Base.Test

@testset begin
    @test PolyaGamma.mean(PolyaGamma.PolyaGammaDistribution(2.0, 3.0)) ≈ 0.3017160845482888
    @test PolyaGamma.mean(PolyaGamma.PolyaGammaDistribution(1.0, 1.0)) ≈ 0.23105857863000487
end
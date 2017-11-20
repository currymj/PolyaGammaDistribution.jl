using PolyaGamma
using Base.Test

@test PolyaGamma.mean(PolyaGamma.PolyaGammaDistribution(2.0, 3.0)) ≈ 0.3017160845482888
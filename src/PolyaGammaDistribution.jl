module PolyaGammaDistribution
using Distributions
using Random: AbstractRNG, GLOBAL_RNG, randexp
using StatsFuns: log1pexp
using SpecialFunctions: lgamma

const TRUNC = 0.64
const cutoff = 1 / TRUNC
const TERMS = 20

"""
    cosh(x) = (1+e⁻²ˣ)/(2e⁻ˣ)
    so logcosh(x) = log((1+e⁻²ˣ)) + x - log(2)
                  = x + log1pexp(-2x) - log(2)
"""
function logcosh(x::Real)
    return x + log1pexp(-2x) - log(2)
end 

"""
A Distribution containing the parameters ``b > 0`` and ``c`` for a Pólya-Gamma
distribution ``PG(b, c)``. Note that while in general ``b`` can be real,
samplers implemented here only work for the integral case.
"""
struct PolyaGamma{T<:Integer, U<:Real} <: ContinuousUnivariateDistribution
    b::T
    c::U
end

function jacobi_logpdf(z, x; ntrunc::Int)
    v = zero(x)
    for n in 0:ntrunc
        v += (iseven(n) ? 1 : -1) * acoef(n, x)
    end
    return logcosh(z) -x*z^2/2 + log(v)
end

"""
    The log coefficients of the infinite sum for the density of PG(b, 0).
    See Polson et al. 2013, section 2.3.
"""
function pg_logcoef(x, b, n)
    lgamma(n+b)-lgamma(n+1) + log(2n+b) -log(2π * x^3)/2 - (2n+b)^2/8x
end
"""
   log density of the PG(b, 0) distribution.
    See Polson et al. 2013, section 2.3.
"""
function pg0_logpdf(x, b; ntrunc::Int)
    v = zero(x)
    for n in 0:ntrunc
        v += (iseven(n) ? 1 : -1) * exp(pg_logcoef(x, b, n))
    end
    return (b-1)*log(2) - lgamma(b) + log(v)
end
"""
    log density of the PG(b, c) distribution.
    See Polson et al. 2013, section 2.2 and equation (5).
"""
function pg_logpdf(b, c, x; ntrunc::Int)
    b*logcosh(c/2) - x*c^2/2 + pg0_logpdf(x, b; ntrunc=ntrunc)
end

function Distributions.logpdf(d::PolyaGamma, x::Real; ntrunc::Int=TERMS)
    if d.b == 1
        return jacobi_logpdf(d.c/2, 4*x; ntrunc=ntrunc) + log(4)
    else
        return pg_logpdf(d.b, d.c, x; ntrunc=ntrunc)
    end
end
Distributions.pdf(d::PolyaGamma, x::Real; ntrunc::Int=TERMS) = exp(Distributions.logpdf(d, x; ntrunc=ntrunc))

"""
Analytically computes the mean of the given PG distribution, using the formula:

``
\frac{b}{2c} \tanh(\frac{c}{2})
``
"""
function Distributions.mean(d::PolyaGamma)
    (d.b / (2.0*d.c)) * tanh(d.c / 2.0)
end

Distributions.rand(d::PolyaGamma) = rand(GLOBAL_RNG, d)
function Distributions.rand(rng::AbstractRNG, d::PolyaGamma)
    rpg_devroye(rng, d.c, d.b, 1)[1]
end

# https://stats.stackexchange.com/questions/122957/what-is-the-variance-of-a-polya-gamma-distribution
# thanks to this nerd for saving me the time of doing the derivation

"""
Analytically computes the variance of the given PG distribution, using the formula

``
\\frac{b}{4c^3} (\\sinh(c) - c) \\sech(\\frac{c}{2})^2
``
"""
function Distributions.var(d::PolyaGamma)
    (d.b / (4 * d.c^3)) * (sinh(d.c) - d.c) * (sech(d.c/2)^2)
end

# functions below are essentially translated from the BayesLogit R package

# cdf of Inverse Gaussian, already helpfully given to us
pigauss(x, μ, λ) = cdf(InverseGaussian(μ, λ), x)

function rtigauss(rng::AbstractRNG, zin, r=TRUNC)
    z = abs(zin)
    μ = 1/z
    x = r + 1
    if (μ > r)
        α = 0.0
        while rand(rng) > α
            ee = randexp(rng, 2)
            while ee[1]^2 > (2 * ee[2] / r)
                ee = randexp(rng, 2)
            end

            x = r / (1 + r * ee[1])^2
            α = exp(-0.5 * z^2 * x)
        end
    else
        while x > r
            λ = 1.0
            y = rand(rng, Normal())^2
            x = μ + 0.5*μ^2 / λ * y - 0.5 * μ / λ * sqrt(4 * μ * λ * y + (μ * y)^2)
            if rand(rng) > (μ/(μ + x))
                x = μ^2/x
            end
        end
    end
    x
end

function mass_texpon(z, x=TRUNC)
    fz = pi^2 / 8 + z^2 / 2
    b = sqrt(1.0 / x) * (x * z - 1)
    a = -1.0 * sqrt(1.0 / x) * (x * z + 1)

    x0 = log(fz) + fz * x
    xb = x0 - z + logcdf(Normal(0,1), b)
    xa = x0 + z + logcdf(Normal(0,1), a)

    qdivp = 4 / pi * (exp(xb) + exp(xa))

    1.0 / (1.0 + qdivp)
end

function rpg_gammasum(num=1, n=1, z=0.0, trunc=200)
    ci = ((1:num).-(1/2)).^2 * pi^2 * 4
    ai = ci + z.^2
    w = zeros(ci)
    for i=1:num
        w[i] = 2.0 * sum(rand(rng, Gamma(n),trunc)./ai)
    end
    w
end

function acoef(n, x, r=TRUNC)
    if ( x>TRUNC )
        pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
    else
        (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
    end
end

# this is the sampler you want for a single element,
# for 1, z
function rpg_devroye_1(rng::AbstractRNG, z::Float64)
    z = abs(z) * 0.5
    fz = pi^2 / 8 + z^2 / 2

    numtrials = 0
    totaliter = 0
    x = 0.0
    while true
        numtrials += 1
        if rand(rng) < mass_texpon(z)
            x = TRUNC + randexp(rng) / fz
        else
            x = rtigauss(rng, z)
        end
        s = acoef(0, x)
        y = rand(rng)*s
        n = 0

        while true
            n += 1
            totaliter += 1
            if n % 2 == 1
                s = s - acoef(n, x)
                if y <= s
                    break
                end
            else
                s = s + acoef(n, x)
                if y > s
                    break
                end
            end
        end
        break
        if y <= s
            break
        end
    end
    0.25 * x
end

function rpg_devroye(rng::AbstractRNG, z=0.0, n=1, num=1)

    x = zeros(num)

    for i=1:num
        x[i] = 0
        for j=1:n
            temp = rpg_devroye_1(rng, z)
            x[i] = x[i] + temp
        end
    end
    x

end

function rpg_alt_1(rng::AbstractRNG, z)
    α = 0.0
    x = 0.0
    while (rand(rng) > α)
        x = rpg_devroye_1(rng, 0.0)
        α = exp(-0.5 * (z * 0.5)^2 * x)
    end
    x
end

function rpg_alt(z, num=1)
    Z = [z for _ in 1:num]
    x = zeros(Z)
    for i=1:num
        x[i] = rpg_alt_1(Z[i])
    end
    x
end

export PolyaGamma
end # module

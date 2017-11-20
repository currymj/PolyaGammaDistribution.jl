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

function mass_texpon(z, r=TRUNC)
    edist = Truncated(Exponential(1), 0.0, r)
    pdf(edist, z)
end

function rpg_gammasum(num=1, n=1, z=0.0, trunc=200)
    ci = ((1:num).-(1/2)).^2 * pi^2 * 4
    ai = ci + z.^2
    w = zeros(ci)
    for i=1:num
        w[i] = 2.0 * sum(rand(Gamma(n),trunc)./ai)
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

function rpg_devroye(z::Float64)
    z = abs(z) * 0.5
    fz = pi^2 / 8 + z^2 / 2

    numtrials = 0
    totaliter = 0
    expd = Exponential(1)
    x = 0.0
    while true
        numtrials += 1
        if rand() < mass_texpon(z)
            x = TRUNC + rand(expd) / fz 
        else
            x = rtigauss(z)
        end
        s = acoef(0, x)
        y = rand()*s
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

function rpg_alt(z)
    α = 0.0
    x = 0.0
    while (rand() > α)
        x = rpg_devroye(0)
        α = exp(-0.5 * (z * 0.5)^2 * x)
    end
    x
end

function rpg_alt(z, num=1)
    Z = [z for _ in 1:num]
    x = zeros(Z)
    for i=1:num
        x[i] = rpg_alt(Z[i])
    end
end

export PolyaGamma, mean
end # module

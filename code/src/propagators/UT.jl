# UT.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Unscented Transform implementation for uncertainty propagation
# Notation from "The Unscented Kalman Filter for Nonlinear Estimation", 
#  Eric A. Wan and Rudolph van der Merwe.
# The methods are typed on L (number of symmetric sigma-vectors) so that the 
# Julia compilar can optimize the method assuming that L stays fixed.
#
# Note that we don't implement a full Kalman filter, merely the unscented transform part,
# as we do not intermediate corrections or anything like that.

struct SigmaVectors{L}
    χ::SMatrix{6,2 * L + 1}
    W::SVector{2 * L + 1}
    W0m::Real
end

"""
    SigmaVectors{L}(dist, α, κ, β)

Generates the sigma-vectors and weights for an UT given the distribution and scaling parameters.
"""
function SigmaVectors{L}(dist::MultivariateDistribution, α, κ, β) where {L}
    χ = @SMatrix zeros((6, 2 * L + 1))
    W = @SVector zeros(2 * L + 1)

    λ = α^2 * (L + κ) - L

    μ = mean(dist)
    χ[:, 1] = μ

    W[1] = λ / (L + λ) + (1 - α^2 + β)

    # The "special" weight for mean propagation
    W0m = λ / (L + λ)

    S = sqrt((L + λ) * cov(dist))

    # Symmetric sigma-vectors
    for i in 2:L
        χ[:, i] = μ + S[:, i]
        χ[:, i+L] = μ - S[:, i]

        W[i] = 1 / (2 * (L + λ))
        W[i+L] = W[i]
    end

    return SigmaVectors{L}(χ, W, W0m)
end

"""
    run_ut{L}(p::ForceModel, dist::MultivariateDistribution, Δt; reltol, abstol, α, κ, β)

Runs the UT propagation for the given distribution, interval of time Δt, with 2L + 1 sigma-points
and α, κ and β scaling parameters (according to "The Unscented Kalman Filter for Nonlinear Estimation").
Note this doesn't implement a full UKF, just the unscented transform part.

Returns the resulting multivariate normal distribution. Note that wathever non-Gaussian dist was used,
the result is always Gaussian!
"""
function run_ut{L}(
    p::ForceModel,
    dist::MultivariateDistribution,
    Δt,
    reltol=1e-10,
    abstol=1e-10,
    α=1e-3,
    κ=0,
    β=2
) where {L}

    sigma = SigmaVectors{L}(dist, α, κ, β)

    # The sigma-vectors are propagated through the non-linear function 
    Threads.@threads for i in 1 .. (2 * L)
        sigma.χ[:, i] = propagate_orbit(p, sigma.χ[:, i], Δt, reltol=reltol, abstol=abstol)
    end

    # Mean is computed, using the "special" weight for the central point
    μ = sigma.χ[:, 1] * sigma.W0m + sum(sigma.χ[:, i] * sigma.W[i] for i in 2:(2*L))

    # Covariance matrix can be computed from the deviation matrix
    dx = sigma.χ .- μ
    P = dx * Diagonal(sigma.W) * dx'

    return MvNormal(μ, P)

end
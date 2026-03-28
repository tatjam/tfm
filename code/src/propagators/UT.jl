# UT.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Unscented Transform implementation for uncertainty propagation
# Notation from "The Unscented Kalman Filter for Nonlinear Estimation", 
#  Eric A. Wan and Rudolph van der Merwe.
#
# Note that we don't implement a full Kalman filter, merely the unscented transform part,
# as we do not intermediate corrections or anything like that.

"""
    SigmaVectors(χ, W, W0m)
"""
struct SigmaVectors{T<:Real}
    # 6-dimensional, 13 sigma points
    χ::SMatrix{6,13,T}
    # One weight for each sigma point
    W::SVector{13,T}
    # Additional weight for average sigma point
    W0c::T
end

"""
    SigmaVectors(μ, P, α, κ, β)

Generates the sigma-vectors and weights for an UT given the distribution mean and covariance matrix, 
and scaling parameters α, κ and β.
"""
function SigmaVectors(μ::AbstractVector{T}, P::AbstractMatrix{T}, α, κ, β) where {T}
    L = 6
    χ = @MMatrix zeros(6, 2 * L + 1)
    W = @MVector zeros(2 * L + 1)

    λ = α^2 * (L + κ) - L

    χ[:, 1] = μ
    W[1] = λ / (L + λ)

    # The "special" weight for mean propagation
    W0c = λ / (L + λ) + (1 - α^2 + β)

    S = sqrt((L + λ) * P)

    # Symmetric L sigma-vectors
    for i in 2:(L+1)
        χ[:, i] = μ + S[:, i-1]
        χ[:, i+L] = μ - S[:, i-1]

        W[i] = 1 / (2 * (L + λ))
        W[i+L] = W[i]
    end

    return SigmaVectors{T}(SMatrix(χ), SVector(W), W0c)
end

"""
    run_ut(p::ForceModel, μ, P, Δt; reltol, abstol, α, κ, β)

Runs the UT propagation for the given distribution, interval of time Δt,
and α, κ and β scaling parameters (according to "The Unscented Kalman Filter for Nonlinear Estimation").
Note this doesn't implement a full UKF, just the unscented transform part.

Returns the resulting multivariate normal distribution.
"""
function run_ut(
    p::ForceModel,
    μ::AbstractVector{T},
    P::AbstractMatrix{T},
    Δt,
    reltol=1e-10,
    abstol=1e-10,
    α=1e-3,
    κ=0,
    β=2
) where {T}

    L = 6
    sigma = SigmaVectors(μ, P, α, κ, β)

    endpoints = Vector{SVector{6,T}}(undef, 2 * L + 1)

    # The sigma-vectors are propagated through the non-linear function 
    Threads.@threads for i in 1:(2*L+1)
        endpoints[i] = propagate_orbit(p, sigma.χ[:, i], Δt, reltol=reltol, abstol=abstol)
    end

    # Mean is computed as simply the mean of all points
    μend = sum(endpoints[i] * sigma.W[i] for i in 1:(2*L+1))

    # Covariance matrix can be computed from the deviation matrix, but we have to use 
    # the special weight for the mean point
    dx = reduce(hcat, endpoints) .- μend
    W = [sigma.W0c; sigma.W[2:end]]

    Pend = dx * Diagonal(W) * dx'

    return MvNormal(μend, Symmetric(Pend))

end
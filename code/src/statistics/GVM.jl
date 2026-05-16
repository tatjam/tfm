# GVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Gauss Von-Mises distribution implementation, a natural representation
# for uncertainty in MEE coordinates. We use notation as defined in
# "Gauss von Mises Distribution for Improved Uncertainty Realism in
#  Space Situational Awareness", Joshua T. Horwood and Aubrey B. Poore, 2014.
# 
# For pure keplerian coordinates (a, e, i, ω, Ω, ν), this distribution is also
# appropriate, but careful choice of i = ω = Ω = 0 has to be done if these
# coordinates are to be considered Gaussian.

struct GaussVonMises{T<:Real,V<:AbstractVector{T},M<:AbstractMatrix{T}}
    # μ mean vector for the Gaussian in Euclidean space
    μ::V
    # α mean vector for the Von Mises (scalar)
    α::T
    # β is a 1-form: ℝⁿ -> ℝ, acting on vectors via dot(β, v)
    # It maps each vector to the effect it has on the angular input to VM
    β::V
    # Γ is a symmetric bilinear form, inducing a quadratic form ℝⁿ -> ℝ, via dot(v, Γ, v)
    Γ::M
    # κ is a scalar that shapes the Von Mises transform, it's non geometric
    κ::T
    # A is a linear map: ℝⁿ -> ℝⁿ, acting on vectors via the matrix product A*v,
    # mapping a "canonical unit sphere" to a "real ellipsoid" for the Gaussian 
    A::LowerTriangular{T}
end

"""
    GaussVonMises(μ, P, α, β, Γ, κ)

Create a Gauss Von Mises distribution for n Euclidean dimensions and a single angular dimension.
The distribution is defined as N(μ, P) × VM(α + β'z + 1/2 z'Γz, κ). A z = x - μ,
where A is the lower triangular Cholesky decomposition of P, and:

     μ is a n dimensional vector representing the mean of the Euclidean elements
     P is a n×n matrix representing the covariance of the Euclidean elements. It's thus symmetric and positive definite.
     α is a real number representing the "angular mean" of the angular element
     β is a n dimensional vector capturing the linear correlation between the euclidean elements and the angular element
     Γ is a n×n matrix used as a quadratic form to tune the shape of the distribution.
     κ is a real number representing the angular concentration of the distribution.

The distribution assigns probabilities to a random tuple (x, θ), where x is the euclidean random variable and θ is the angular random variable.
    
"""
function GaussVonMises(μ::V, P::M, α::T, β::V, Γ::M, κ::T) where {T<:Real,V<:AbstractVector{T},M<:AbstractMatrix{T}}
    chol = cholesky(Symmetric(P)).L
    GaussVonMises(μ, α, β, Γ, κ, chol)
end

"""
    decanonicalize(dist, v)

Transforms a random vector distributed under a canonical GVM (i.e. GVM(0, I, 0, 0, 0, κ)) to its corresponding non-canonical GVM dist, with same κ.

The vector is of the form [euclidean part..., angular value]!
"""
function decanonicalize(dist::GaussVonMises{T, V, M}, v::AbstractVector{T}) where {T, V, M}

    ceuc = v[1:(length(v) - 1)]
    cang = v[length(v)]

    euc = dist.μ + dist.A * ceuc
    ang = cang + dist.α + dot(dist.β, ceuc) + 0.5 * dot(ceuc, dist.Γ, ceuc)

    SA[euc..., ang]
end

"""
    rand(rng, d::GaussVonMises)

Sample a GaussVonMises distribution, returning the (x, θ) sample
"""
function Base.rand(rng::AbstractRNG, d::GaussVonMises)
    # z is distributed as N(0, Id), thus is in our "canonical" space
    z = randn(rng, eltype(d.μ), length(d.μ))
    # x is distributed as N(μ, P), using A to go from canonical to real coordinates
    x = d.μ + d.A * z

    # θ is constructed from α, the action of β on z, and the action of Γ in z
    Θ = d.α + dot(z, d.β) + 0.5 * dot(z, d.Γ, z)

    # Note: Distributions.jl VonMises is always centered on the mean, so
    # we "unwrap" so it lives on [-π, π)
    vm = mod(rand(rng, VonMises(Θ, d.κ)) + π, 2 * π) - π
    (x, vm)
end

function Base.rand(rng::AbstractRNG, d::GaussVonMises, dims::NTuple{N, Int}) where N
    reshape([rand(rng, d) for _ in 1:prod(dims)], dims)
end

function mahalanobis(x::V, θ::T, dist::GaussVonMises{T, V, M})
    where {T<:Real,V<:AbstractVector{T},M<:AbstractMatrix{T}}

    deuclid = x - dist.μ
    z = dist.A \ deuclid 

    # Alternative way to cheaply compute (x-μ)ᵀ(AAᵀ)⁻¹(x-μ)
    # z = A⁻¹ (x-μ), thus
    # zᵀz = (x-μ)ᵀ A⁻ᵀ A⁻¹ (x-μ) = (x-μ)ᵀ(AAᵀ)⁻¹(x-μ)
    euclid_term = dot(z, z)

    expected_ang =  dist.α + dot(dist.β, z) + 0.5 * dot(z, dist.Γ, z)
    dang = θ - expected_ang
    t1 = sin(0.5 * dang)
    ang_term = T(4.0) * dist.κ * t1 * t1

    return euclid_term + ang_term
end

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
    nearest_pd_matrix(mat)

Obtains the nearest positive definite matrix to mat, by clamping the eigen values.

WARNING: Badly behaved if mat has both large and very small eigenvalues.
"""
function nearest_pd_matrix(mat)
    λ, V = eigen(Symmetric(mat))
    # minλ = 1e-6 * maximum(λ)
    minλ = eps(eltype(mat)) * maximum(λ)

    for λi in λ
        if λi < minλ
            @warn "Very small eigen value $(λi) clamped to $(minλ)"
        end
    end
    return Symmetric(V * Diagonal(max.(λ, minλ)) * V')
end

"""
    equilibrated_sqrt(P; λmin = nothing)

Return D and r, with `S = Dr => S * S' ≈ P`, using a scaled space defined as a scaling of the
coordinate axes by their standard deviations (if P is a covariance matrix) to improve
numerical precision.

If λmin is not given, uses a small numerically meaningful minimum value, chosen to be
a bit above the numerical tolerance of the eigen value decomposition precision.

"""
function equilibrated_sqrt(P::AbstractMatrix; λmin = nothing)
    s = sqrt.(max.(diag(P), 0.0))
    # If any element is zero or negative, keep it unscaled
    s = map(x -> x > 0 ? x : one(x), s)

    D = Diagonal(s)
    # The left division by D is scaling of the rows, and the right division of the
    # columns. Thus, C = D⁻¹ P D⁻¹ and P = D C D
    C = Symmetric(D \ Matrix(P) / D)
    F = eigen(C)

    tol = something(λmin, 10 * size(C, 1) * eps(float(eltype(C))) * maximum(F.values))

    λ = max.(F.values, tol)

    # Then C = r * rᵀ 
    r = F.vectors * Diagonal(sqrt.(λ))
    # thus (D * r) (D * r)ᵀ = D C D = P
    return D, r
end

"""
    equilibrated_pdm(raw; λmin = nothing)

Builds a `PDMat` from a covariance matrix `raw`, by clamping its negative eigen-values
in a scaled space to guarantee better numerical behavior.

If λmin is not given, uses a small numerically meaningful minimum value, chosen to be
a bit above the numerical tolerance of the eigen value decomposition precision.
"""
function equilibrated_pdm(raw::AbstractMatrix; λmin = nothing)
    D, r = equilibrated_sqrt(raw; λmin = λmin)
    # The restored P matrix cholesky factor, computing cholesky of the scaled matrix
    L_full = D * Matrix(cholesky(r * r').L)
    # By definition of Cholesky factoring
    P_restore = L_full * L_full'
    return PDMat(Symmetric(L_full * L_full'), Cholesky(L_full, 'L', 0))
end

"""
    SigmaVectors(χ, W, W0c)
"""
struct SigmaVectors{L,N,T<:Real}
    # L-dimensional, 2L+1 sigma points
    χ::SMatrix{L,N,T}
    # One weight for each sigma point
    W::SVector{N,T}
    # Additional weight for average sigma point
    W0c::T
end

"""
    SigmaVectors(μ, P, α, κ, β; equilibrate = false)

Generates the sigma-vectors and weights for an UT given the distribution mean
and covariance matrix, and scaling parameters α, κ and β.

With `equilibrate = true` the square root of `P` is computed in an equilibrated
space (see `equilibrated_sqrt`).
"""
function SigmaVectors(
    μ::AbstractVector{T},
    P::AbstractMatrix{T},
    α,
    κ,
    β;
    equilibrate = false,
) where {T}
    L = length(μ)
    N = 2 * L + 1
    χ = @MMatrix zeros(L, N)
    W = @MVector zeros(N)
    λ = α^2 * (L + κ) - L

    χ[:, 1] = μ
    W[1] = λ / (L + λ)

    # The "special" weight for mean propagation
    W0c = λ / (L + λ) + (1 - α^2 + β)

    S = if equilibrate
        first(equilibrated_sqrt((L + λ) * Matrix(P)))
    else
        sqrt((L + λ) * P)
    end

    # Symmetric L sigma-vectors
    for i = 2:(L+1)
        χ[:, i] = μ + S[:, i-1]
        χ[:, i+L] = μ - S[:, i-1]
        W[i] = 1 / (2 * (L + λ))
        W[i+L] = W[i]
    end
    return SigmaVectors{L,N,T}(SMatrix(χ), SVector(W), W0c)
end

"""
    ut_propagate(f, μ, P; α, κ, β, equilibrate = false) 

Propagates a normal distribution through function f, and reconstruct the normal distribution
at the other side by means of the Unscented Transform.

With `equilibrate = true` both the sigma-point generation and the output
covariance repair are done in an equilibrated space (see `equilibrated_sqrt`
and `equilibrated_pdm`).
"""
function ut_propagate(
    f,
    μ::AbstractVector{T},
    P::AbstractMatrix{T};
    α = 1e-3,
    κ = 0,
    β = 2,
    equilibrate = false,
) where {T}
    L = length(μ)
    N = 2 * L + 1
    sigma = SigmaVectors(μ, P, α, κ, β; equilibrate = equilibrate)
    endpoints = Vector{SVector{L,T}}(undef, N)
    # The sigma-vectors are propagated through the non-linear function 
    Threads.@threads for i = 1:N
        endpoints[i] = f(sigma.χ[:, i])
    end

    # Mean is computed as simply the mean of all points
    μend = sum(endpoints[i] * sigma.W[i] for i = 1:N)
    # Covariance matrix can be computed from the deviation matrix, but we have to use 
    # the special weight for the mean point
    dx = reduce(hcat, endpoints) .- μend
    W = [sigma.W0c; sigma.W[2:end]]

    Pend = if equilibrate
        equilibrated_pdm(dx * Diagonal(W) * dx')
    else
        nearest_pd_matrix(dx * Diagonal(W) * dx')
    end

    return MvNormal(μend, Pend)
end

"""
    run_ut(p::ForceModel, μ, P, Δt; reltol, abstol, α, κ, β, equilibrate)

Runs the UT propagation for the given distribution, interval of time Δt,
and α, κ and β scaling parameters (according to "The Unscented Kalman Filter for Nonlinear Estimation").
Note this doesn't implement a full UKF, just the unscented transform part.

With `equilibrate = true` the covariance is equilibrated before any
eigendecomposition or Cholesky (see `ut_propagate`), which is required for
covariances mixing units of very different magnitude (e.g. MEE elements).

Returns the resulting multivariate normal distribution.
"""
function run_ut(
    p::ForceModel,
    μ::AbstractVector{T},
    P::AbstractMatrix{T},
    Δt,
    reltol = 1e-10,
    abstol = 1e-10,
    α = 1e-3,
    κ = 0,
    β = 2;
    equilibrate = false,
) where {T}
    return ut_propagate(
        v -> propagate_orbit(p, v, Δt, reltol = reltol, abstol = abstol),
        μ,
        P,
        α = α,
        κ = κ,
        β = β,
        equilibrate = equilibrate,
    )
end

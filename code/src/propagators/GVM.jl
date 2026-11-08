# GVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Gauss Von-Mises propagator, essentially the unscented transform adapted to GVM distribution.
# Implemented according to 
# "Gauss von Mises Distribution for Improved Uncertainty Realism in
#  Space Situational Awareness", Joshua T. Horwood and Aubrey B. Poore, 2014.


"""
L-dimensional, 2L+1 sigma points. If n is the number of Euclidean dimensions, 
GVM lives on ℝⁿ ⨯ S, so L = n + 1 and N = 2L + 1 = 2n + 3 sigma points.

Indices are:
    0      -> Central point
    1:2    -> Angular offset
    3:2L+1 -> Euclidean offsets
"""
struct GVMSigmaVectors{L,N,T<:Real}
    χ::SMatrix{L,N,T}  # Sigma points matrix
    W::SVector{N,T}     # Weights for each sigma point
end

"""
   b12(κ) = (1 - I₁ / I₀, 1 - I₂ / I₀)

Evaluates modified Bessel functions of the first kind at κ.
"""
function b12(κ::T) where {T}
    # Note, because we always compute their division, we can use
    # the scaled bessel functions which don't explode with κ
    bessels = map(n -> Bessels.besselix(n, κ), 0:2)
    return (T(1) - bessels[2] / bessels[1], T(1) - bessels[3] / bessels[1])
end

"""
    sigma_params(κ)

Returns acos(b₂(κ) / (2 b₁(κ)) - 1) and b₁²(κ) / (4 b₁(κ) - b₂(κ)) using the asymptotic approximation for them when κ is large.

See annex for derivation of the asymptotic expressions.
"""
function sigma_params(κ::T) where {T}
    if κ > 1e5
        return sqrt(3 / κ), T(1) / 6
    else
        b1, b2 = b12(κ)
        return acos(b2 / (T(2) * b1) - T(1)), (b1 * b1) / (T(4) * b1 - b2)
    end
end

function GVMSigmaVectors(dist::GaussVonMises{T}) where {T}
    n = length(dist.μ)
    L = n + 1
    N = 2L + 1

    # Canonical sigma vectors depend solely on κ
    ξc = sqrt(T(3))
    ηc, wcη = sigma_params(dist.κ)
    wcξ = T(1) / T(6)
    wc0 = T(1) - T(2) * wcη - T(2) * n * wcξ

    # Generate canonical sigma vectors
    χ = begin
        χbuf = MMatrix{L,N,T}(undef)
        χbuf[:, 1] = SA[zero(dist.μ)..., T(0)]   # χ00
        χbuf[:, 2] = SA[zero(dist.μ)..., T(ηc)] # χη0
        χbuf[:, 3] = SA[zero(dist.μ)..., -T(ηc)] # χη1

        # Symmetric vectors around the origin for the Euclidean part
        for i = 1:n
            idx = 3 + i
            χbuf[:, idx] = SA[zero(dist.μ)..., T(0)]
            χbuf[i, idx] = ξc
            χbuf[:, idx+n] = SA[zero(dist.μ)..., T(0)]
            χbuf[i, idx+n] = -ξc
        end

        reduce(hcat, decanonicalize.(Ref(dist), eachcol(SMatrix(χbuf))))
    end

    w = SA[wc0, wcη, wcη, fill(wcξ, 2n)...]
    return GVMSigmaVectors(SMatrix(χ), SVector(w))
end

"""
    GVMLeastSquares{S, L, N, M, T}

Residual struct to minimize via least squares on the GVM angular parameters (α, β, Γ),
the difference of the Mahalanobis-Von Mises of the sigma-points to the initial
distribution (lₑ), and their distance after the transformation (lₐ).

Note that NLLS always minimizes the absolute value!

Type parameters:
    - S: Number of free parameters in Γ, S = ½L(L-1)
    - L: Problem dimensions, L = n + 1
    - N: Number of sigma points, N = 2L + 1
    - M: Number of euclidean dimensions (n)
"""
struct GVMLeastSquares{S,L,N,M,T} <: NLLSsolver.AbstractResidual
    # Initial mahalanobis distance of each sigma-point
    lₑ::SVector{N,T}
    end_σ::SMatrix{L,N,T}
    κ::T
    end_μ::SVector{M,T}
    end_A::LowerTriangular{T,SMatrix{M,M,T}}
end

# NLLSsolver Boilerplate
Base.eltype(::GVMLeastSquares{S,L,N,M,T}) where {S,L,N,M,T} = T
# How many variables do we optimize? (α, β, Γ)
NLLSsolver.ndeps(::GVMLeastSquares) = static(3)
# How many residuals do we minimize?
NLLSsolver.nres(::GVMLeastSquares{S,L,N,M,T}) where {S,L,N,M,T} = static(N)
# See below, (α, β, Γ) -> indices of addvariable from top to bottom
NLLSsolver.varindices(::GVMLeastSquares) = SVector(1, 2, 3)

function NLLSsolver.getvars(::GVMLeastSquares{S,L,N,M,T}, vars::Vector) where {S,L,N,M,T}
    return (
        vars[1]::NLLSsolver.EuclideanVector{1,T},  # end_α
        vars[2]::NLLSsolver.EuclideanVector{M,T},  # end_β
        vars[3]::NLLSsolver.EuclideanVector{S,T},  # end_Γ (flattened)
    )
end

function NLLSsolver.computeresidual(
    res::GVMLeastSquares{S,L,N,M,T},
    end_α,
    end_β,
    end_Γ_vec,
) where {S,L,N,M,T}
    # NLLS uses dual numbers, so extract the type if any of the variables
    # is being differentiated on
    T_AD = promote_type(eltype(end_α), eltype(end_β), eltype(end_Γ_vec))
    n_Γ = L - 1

    end_Γ = begin
        tmp = zeros(T_AD, n_Γ, n_Γ)
        tmp[triu(ones(Bool, n_Γ, n_Γ))] = end_Γ_vec
        Symmetric(tmp)
    end

    # Constants as dual types so we don't get typing errors / zero derivatives
    μ_ad = T_AD.(res.end_μ)
    A_ad = LowerTriangular(T_AD.(res.end_A))
    κ_ad = T_AD(res.κ)

    # α is stored as a vector of one entry
    α_scalar = end_α isa AbstractArray ? end_α[1] : end_α

    end_dist = GaussVonMises(μ_ad, α_scalar, end_β, end_Γ, κ_ad, A = A_ad)

    return SVector{N,T_AD}(ntuple(Val(N)) do i
        # Original vs ending mahalanobis, |lₑ - lₐ| is minimized.
        T_AD(res.lₑ[i]) - mahalanobis(res.end_σ[:, i], end_dist)
    end)
end

"""
    GVMRegularizer{S, M, T}

We have 1 {α} + L - 1 {β} + ½L(L-1) {Γ} = ½L(L+1) variables to optimize,
and 2L+1 equations (sigma-points). Thus, there are many valid solutions as
long as L >= 4 (for L = 6, we have 8 more unknowns than equations). To try
to keep the solution close to the initial estimate, we add this regularizer
which just gives cost to diverging from the initial estimate via a weight.

The weight is chosen to be √(1+κ), this is so the weight of the regularizer
is more or less "equally important" as the mahalanobis optimization. This value
comes from the GVM mahalanobis distance formula.

First, consider that the sigma-points are chosen precisely to lay at a constant
mahalanobis distance from the mean. Consider that each sigma-point is at a
angular distance "of order φ" from the mean, the mahalanobis distance's
angular term is thus 4κ sin²(0.5 φ), which for small φ yields κφ².

Then, φ = O(1/√κ), which we assume to hold after the nonlinear propagation.

Now, consider the derivative of the residual of a sigma-point from the solver
wrt. its angular part:
    ∂l/∂φ = ∂/∂φ (lₑ - zᵀz - 4κ sin²(0.5 φ) ) = - 2 κ sin(φ)

The Jacobian "J" has columns for one sigma-point:
    1 column for        ∂l/∂α = ∂l/∂φ ∂φ/∂α
    L - 1 columns for   ∂l/∂β = ∂l/∂φ ∂φ/∂β
    ½L(L-1) columns for ∂l/∂Γ = ∂l/∂φ ∂φ/∂Γ

Where we see the ∂l/∂φ = 2κ sin(φ) is scaling every term in the
Jacobian matrix. For small angles, we have a scaling O(κ / √κ) = O(√κ).

Finally, for the regularizer to remain O(1) at small κ, we add the plus one,
giving w = √(1+κ).
"""
struct GVMRegularizer{S,M,T} <: NLLSsolver.AbstractResidual
    α₀::T
    β₀::SVector{M,T}
    Γ₀::SVector{S,T}
    w::T
end

Base.eltype(::GVMRegularizer{S,M,T}) where {S,M,T} = T
NLLSsolver.ndeps(::GVMRegularizer) = static(3)
NLLSsolver.nres(::GVMRegularizer{S,M,T}) where {S,M,T} = static(1 + M + S)
NLLSsolver.varindices(::GVMRegularizer) = SVector(1, 2, 3)

function NLLSsolver.getvars(::GVMRegularizer{S,M,T}, vars::Vector) where {S,M,T}
    return (
        vars[1]::NLLSsolver.EuclideanVector{1,T},
        vars[2]::NLLSsolver.EuclideanVector{M,T},
        vars[3]::NLLSsolver.EuclideanVector{S,T},
    )
end

function NLLSsolver.computeresidual(reg::GVMRegularizer{S,M,T}, α, β, Γ_vec) where {S,M,T}
    α_scalar = α isa AbstractArray ? α[1] : α
    return SVector(
        reg.w * (α_scalar - reg.α₀),
        (reg.w .* (β .- reg.β₀))...,
        (reg.w .* (Γ_vec .- reg.Γ₀))...,
    )
end

function gvm_propagate(
    f,
    dist::GaussVonMises{T};
    slopes::Symbol = :ad,
    refine::Bool = true,
) where {T}
    n = length(dist.μ)
    L = n + 1
    N = 2L + 1
    S = (L * (L - 1)) ÷ 2

    # Step 1: Generate and propagate sigma vectors
    sigma = GVMSigmaVectors(dist)
    endpoints = Vector{SVector{L,T}}(undef, N)

    Threads.@threads for i = 1:N
        endpoints[i] = f(sigma.χ[:, i])
    end

    # Step 2: GVM quadrature for the Euclidean part.
    end_μ, end_P, end_A = let
        out_μ = sum(endpoints .* sigma.W)[1:n]

        dx = reduce(hcat, endpoints)[1:n, :] .- out_μ

        out_P = equilibrated_pdm(dx * Diagonal(sigma.W) * dx')
        out_A = cholesky(Symmetric(out_P)).L
        out_μ, out_P, out_A
    end

    # Step 3: Initial estimate of α, β, Γ and κ by composing the prior angular
    # mean θ(z) = α + βᵀz + ½zᵀΓz (evaluated at the sigma points) with f, to
    # second order in the prior canonical coordinate z, and inverting the
    # Euclidean output map to express the result in the end canonical frame.
    #
    # The slopes come either from AD at the central point (:ad), or from the
    # propagated nodes themselves (:nodes): statistical linear regression, the
    # UKF cross-covariance principle via Stein's lemma. The ±ξeᵢ node pairs
    # give slopes and diagonal curvatures that are exact for quadratic maps,
    # and the ±η pair gives the φ-slope through the von Mises identity
    # E[∂φF] = κE[sinφ·F]. The node version needs no derivative passes through
    # the ODE, but cannot recover off-diagonal curvature.
    end_α, end_β, end_Γ, end_κ = let
        grad, Q, M, q, c₀, g = if slopes === :ad
            # Proceed via automatic differentation
            central_point = sigma.χ[:, 1]

            δf = ForwardDiff.jacobian(f, central_point)
            # δ²f[i, j, k] = ∂²f_i/∂x_j∂x_k
            δ²f = reshape(
                ForwardDiff.jacobian(x -> vec(ForwardDiff.jacobian(f, x)), central_point),
                L,
                L,
                L,
            )

            # Prior canonical Jacobian of the state: s(z) ≈ c + Jz
            J = vcat(dist.A, permutedims(dist.β))

            grad = J' * δf[L, :]
            Q = δf[L, L] * dist.Γ + J' * δ²f[L, :, :] * J
            M = δf[1:n, :] * J
            q = ntuple(i -> J' * δ²f[i, :, :] * J + δf[i, L] * dist.Γ, Val(n))
            (grad, Q, M, q, δf[L, L], δf[1:n, L])
        elseif slopes === :nodes
            # Proceed via a similar process to the UKF (see writeup)
            ξ = sqrt(T(3))
            η, wη = sigma_params(dist.κ)
            F0 = endpoints[1]
            Dz = reduce(hcat, ((endpoints[3+i] - endpoints[3+n+i]) / (2ξ) for i = 1:n))
            D2 = reduce(hcat, ((endpoints[3+i] + endpoints[3+n+i] - 2F0) / ξ^2 for i = 1:n))
            Dφ = dist.κ * wη * sin(η) * (endpoints[2] - endpoints[3])
            (
                Dz[L, :],
                Diagonal(D2[L, :]),
                Dz[1:n, :],
                ntuple(k -> Diagonal(D2[k, :]), Val(n)),
                Dφ[L],
                Dφ[1:n],
            )
        else
            error("slopes must be :ad or :nodes, got $slopes")
        end

        # Invert the Euclidean map, z = M⁻¹ẋ − ½M⁻¹q(M⁻¹ẋ) with ẋ = x' − f_x(c),
        # and reexpress in the end canonical frame z' = A'⁻¹(x' − out_μ):
        #   Γ' = Pᵀ Qeff P
        #   β' = Pᵀ (grad − Qeff p)
        #   α' = f_θ(c) − β'ᵀz'_c − ½z'_cᵀΓ'z'_c
        # with P = M⁻¹A', p = M⁻¹(f_x(c) − out_μ), h = M⁻ᵀgrad and
        # Qeff = Q − Σᵢ hᵢqᵢ.
        Minv = inv(M)
        h = Minv' * grad
        Qeff = Q - sum(h[i] * q[i] for i = 1:n)

        P = Minv * end_A
        δ = endpoints[1][1:n] - end_μ
        p = Minv * δ

        Γ_est = P' * Qeff * P
        β_est = P' * (grad - Qeff * p)

        # α is referenced to the new canonical origin out_μ: the propagated
        # central point sits at z'_c ≠ 0 in general, and the angular mean
        # surface must pass through it, θ'(f_x(c)) = f_θ(c)
        z_c = end_A \ δ
        α_est = endpoints[1][L] - dot(β_est, z_c) - 0.5 * dot(z_c, Γ_est, z_c)

        # The angular noise ν ~ VM(κ) reaches θ' scaled by
        #   c = ∂f_θ/∂θ − β'ᵀA'⁻¹(∂f_x/∂θ),
        # the direct scaling minus the part the conditional mean surface
        # absorbs through the ν-induced Euclidean motion. Hence κ' = κ/c²
        # (inert when c ≈ 1, e.g. short propagation steps).
        c = c₀ - dot(β_est, end_A \ g)
        end_κ = isfinite(c) && abs(c) > 1e-3 ? dist.κ / c^2 : dist.κ

        α_est, β_est, Γ_est, end_κ
    end

    # Step 4 (optional): refine (α, β, Γ) by least squares
    if refine
        end_α, end_β, end_Γ = let
            lₑ = [mahalanobis(s, dist) for s in eachcol(sigma.χ)]

            problem = NLLSsolver.NLLSProblem()

            Γ_vec = end_Γ[triu(ones(Bool, n, n))]
            NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(end_α))
            NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(end_β...))
            NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(Γ_vec...))

            NLLSsolver.addcost!(
                problem,
                GVMLeastSquares{S,L,N,n,T}(
                    SVector(lₑ),
                    SMatrix{L,N,T}(reduce(hcat, endpoints)),
                    end_κ,
                    SVector{n}(end_μ),
                    LowerTriangular(SMatrix{n,n}(end_A)),
                ),
            )
            NLLSsolver.addcost!(
                problem,
                GVMRegularizer{S,n,T}(
                    end_α,
                    SVector{n}(end_β),
                    SVector{S}(Γ_vec),
                    sqrt(1 + dist.κ),
                ),
            )

            NLLSsolver.optimize!(problem, NLLSsolver.NLLSOptions())

            α_final = problem.variables[1][1]
            β_final = problem.variables[2]
            Γ_final_vec = problem.variables[3]

            Γ_final = begin
                tmp = zeros(T, n, n)
                tmp[triu(ones(Bool, n, n))] = Γ_final_vec
                Symmetric(tmp)
            end

            α_final, β_final, Γ_final
        end
    end

    return GaussVonMises(end_μ, end_α, end_β, end_Γ, end_κ, A = end_A)
end

"""
    run_gvm(p::ForceModel, dist::GaussVonMises, Δt; reltol, abstol, slopes, refine)


Slopes is either :ad (automatic differentaition) or :nodes (via UT like method).
"""
function run_gvm(
    p::ForceModel,
    dist::GaussVonMises{T,V},
    Δt,
    reltol = 1e-10,
    abstol = 1e-10;
    slopes::Symbol = :ad,
    refine::Bool = true,
) where {T,V}
    return gvm_propagate(dist; slopes, refine) do v
        propagate_orbit(p, v, Δt, reltol = reltol, abstol = abstol)
    end
end

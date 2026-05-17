# GVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Gauss Von-Mises propagator, essentially the unscented transform adapted to GVM distribution.
# Implemented according to 
using NLLSsolver: NLLSInternal
# "Gauss von Mises Distribution for Improved Uncertainty Realism in
#  Space Situational Awareness", Joshua T. Horwood and Aubrey B. Poore, 2014.
# 


"""
L-dimensional, 2L+1 sigma points, note that if n is the number of
Euclidean dimensions, GVM lives on ℝⁿ ⨯ S, thus we have L = n + 1
actual variables, and the UT uses 2L+1 = 2n + 3 sigma points.

Note that in this structure L represents the total number of
dimensions (n + 1), and N = 2L + 1, we don't use the euclidean dimensions n.
"""
struct GVMSigmaVectors{L,N,T<:Real}
    # χ[:, 1] is the central point "00",
    # χ[:, 2] is the positive angular offset point "η0"
    # χ[:, 3] is the negative angular offset point "η1"
    # χ[:, 4..N] are the 2n symmetric euclidean offset points "0ξ"
    χ::SMatrix{L,N,T}
    # Weights for each sigma point
    W::SVector{N,T}
end

"""
   b12(κ) = (1 - I₁ / I₀, 1 - I₂ / I₀)

Where I₂, I₁ and I₀ are the modified Bessel functions of the first kind evaluated at κ.

Returns the pair as it's more computationally efficient to compute all Bessel functions at once.
"""
function b12(κ::T) where {T}
    bessels = Bessels.besseli(0:2, κ)
    return (1 - bessels[2] / bessels[1], 1 - bessels[3] / bessels[1])
end

"""
   GVMSigmaVectors(dist::GaussVonMises) 

Generates the sigma-vectors and weights for a GVM "UT" given the starting distribution.

Proceeds via the canonical GVM (i.e. normalized distribution) as described in the paper.
"""
function GVMSigmaVectors(dist::GaussVonMises{T}) where {T}
    n = length(dist.μ)
    L = n + 1
    N = 2 * L + 1

    # The canonical sigma vectors depend only on the value of κ of the original
    # distributions
    (b1, b2) = b12(dist.κ)

    ξc = sqrt(T(3))
    ηc = acos(b2 / (T(2) * b1) - T(1))

    wcξ = T(1) / T(6)
    wcη = (b1 * b1) / (T(4) * b1 - b2)
    wc0 = T(1) - T(2) * wcη - T(2) * n * wcξ

    # Generate the canonical sigma vectors, living in the canonical
    # GVM with κ = dist.κ
    χ = begin
        χbuf = MMatrix{L,N,T}(undef)
        χbuf[:, 1] = SA[zero(dist.μ)..., T(0)]  # χ00
        χbuf[:, 2] = SA[zero(dist.μ)..., T(ηc)] # χη0
        # Small typo on paper for the next one (there are two
        # sigma points for the angular direction, not just one)
        χbuf[:, 3] = SA[zero(dist.μ)..., -T(ηc)] # χη1

        # Symmetric vectors around the origin for euclidean part
        for i in 1:n
            idx = 3 + i
            # Positive member
            χbuf[:, idx] = SA[zero(dist.μ)..., T(0)]
            χbuf[i, idx] = ξc
            # Negative member
            χbuf[:, idx+n] = SA[zero(dist.μ)..., T(0)]
            χbuf[i, idx+n] = -ξc
        end

        decanonicalize.(dist, SMatrix(χbuf))
    end

    w = SA[wc0, wcη, wcη, (ones(n) * wcξ)...]

    GVMSigmaVectors(χ, w)
end

"""
GVMLeastSquares{L, N, S, T}(lₑ, end_σ, κ)

Construct a GVMLeastSquares problem given the Mahalanobis distance of each starting sigma-point,
the ending sigma-points end_σ (in ambient space) and κ, assumed to be conserved.

Type parameters:
    L = Dimensionality of the sigma points (n + 1)
    N = Number of sigma points, N = 2L + 1
    S = Number of free parameters in Γ, ½L(L - 1) elements:
        Γ is a symmetrirc n×n matrix, thus it has ½n(n+1), in terms of L
        we can write n = (L - 1), thus it has ½((L-1)² + (L-1)) = ½(L² + 1 - 2L + L - 1) = ½(L² - L) = ½L(L - 1)

Let (xᵢ,θᵢ) be our starting sigma points, and (end_xᵢ,end_θᵢ) our transformed sigma points.
The idea is to find end_α, end_β, end_Γ such that
    ∑ᵢ(lₑ(xᵢ,θᵢ; μ, P, α, β, Γ, κ) - lₐ(end_xᵢ,end_θᵢ; end_μ, end_P, end_α, end_β, end_Γ, end_κ))²
is minimum (assuming end_κ = κ).
lₑ and lₐ are the Mahalanobis distances of the starting and ending points respectively.

Intuitively, we want to find the GVM distribution that's most similar to the original in terms of
the "likelyhood" of the sigma-points being samples of said distribution.
"""
struct GVMLeastSquares{L, N, S, T} <: NLLSsolver.AbstractResidual
    # lₑ, starting Mahalanobis distances of the sigma-points, precomputed
    lₑ::SVector{N, T}
    # (end_xᵢ,end_θᵢ), N ending sigma-points in ambient space
    end_σ::SMatrix{L, N, T}
    # Assumed equal before and after transformation
    κ::T
    # Ending μ from quadrature
    end_μ::AbstractVector{T}
    # Ending A from quadrature
    end_A::LowerTriangular{T}
end

# Boiler plate for NLLSsolver...
Base.eltype(::GVMLeastSquares{L, N, S, T}) where {L, N, S, T} = T
# Number of variables we optimize, 3 as we have {end_α, end_β, end_Γ}
NLLSsolver.ndeps(::GVMLeastSquares{L, N, S, T}) where {L, N, S, T} = static(3)
# We have a single residual 
NLLSsolver.nres(::GVMLeastSquares{L, N, S, T}) where {L, N, S, T} = static(1)
# We store each of the variables at indices 1 2 and 3 
NLLSsolver.varindices(::GVMLeastSquares{L, N, S, T}) where {L, N, S, T} = SVector(1, 2, 3)
# Fetch the variables
NLLSsolver.getvars(::GVMLeastSquares{L, N, S, T}, vars::Vector) where {L, N, S, T} = (
    vars[1]::NLLSsolver.EuclideanVector{1, T}, # end_α
    vars[2]::NLLSsolver.EuclideanVector{N, T}, # end_β
    vars[3]::NLLSsolver.EuclideanVector{S, T}, # end_Γ, stored as the independent elements
)

# Residual is the cost function which will be minimized
function NLLSsolver.computeresidual(res::GVMLeastSquares{L, N, S, T}, end_α, end_β, end_Γ_vec) where {L, N, S, T}
    n_Γ = (L-1)
    end_Γ = let
        tmp = zeros(n_Γ, n_Γ)
        tmp[axes(tmp, 1) .>= axes(tmp, 2)'] = end_Γ_vec
        Symmetric(tmp)
    end

    end_dist = GaussVonMises(res.end_μ, end_α, end_β, end_Γ, res.κ, A=res.end_A)
    sum(zip(lₑᵢ, eachcol(end_σᵢ))) do
        lₐᵢ = mahalanobis(end_σᵢ[1:(L-1)], end_σᵢ[L], end_dist)
        r = lₑᵢ + lₐᵢ

        return r * r
    end
end

function gvm_propagate(
    f,
    dist::GaussVonMises{T,V}
) where {T,V}

    n = length(dist.μ)
    L = n + 1
    N = 2 * L + 1

    # Step 1: generate the sigma vectors and propagate them
    sigma, endpoints = let
        sigma = GVMSigmaVectors(dist)
        endpoints = Vector{SVector{L,T}}(undef, N)

        Threads.@threads for i in 1:N
            endpoints[i] = f(sigma.χ[:, i])
        end

        sigma, endpoints
    end

    # Step 2: GVM cuadrature on the resulting sigma-vectors for the
    #         euclidean part of the distribution
    end_P, end_A = let
        end_μ = sum(endpoints .* sigma.W)
        dx = reduce(hcat, endpoints) .- end_μ
        end_P = nearest_pd_matrix(dx * Diagonal(sigma.W) * dx')
        end_A = cholesky(Symmetric(end_P)).L
        end_P, end_A
    end

    # Step 3: Estimate hα, hβ, hΓ using derivatives of f
    # RESEARCH: It would be wise to experiment using the derivatives that
    # can be inferred from the endpoints instead of a full forward diff of f.
    end_α, end_β, end_Γ = let

        # δf[i, j] = δfᵢ/δxⱼ, so we have...
        δf = ForwardDiff.jacobian(f, sigma.N[:, 1])
        δₓfx = δf[1:(L-1), 1:(L-1)]
        δₓfα = δf[L, 1:(L-1)]

        # δ²fθ[i, j] = δ²f / δxᵢδxⱼ
        δ²fα = ForwardDiff.hessian(x -> f(x)[L], sigma.N[:, 1])
        δ²ₓfα = δ²fα[1:(L-1), 1:(L-1)]

        # Written so they act on canonical vectors
        canon_δₓfx = inv(end_A) * δₓfx * dist.A
        canon_δₓfx = inv(end_A) * δₓfx * dist.A
        canon_δ²ₓfα = dist.A' * δ²ₓfα  * dist.A 

        Δβ = dist.A' * δₓfα
        end_β_original_coords = (dist.β + Δβ)


        # Initial estimates of α, β and Γ
        end_α = endpoints[1][L]
        end_β = inv(end_A) * δₓfx * dist.A * end_β_original_coords
        end_Γ = canon_δₓfx * (dist.Γ + canon_δ²ₓfα) * canon_δₓfx'

        end_α, end_β, end_Γ
    end


    # Step 4: Solve least square problem to refine end_α, end_β, end_Γ.
    # We do so in the "log space", as min(log(x)) = min(x) by monotonicity of log.
    # RESEARCH: Any clever way to avoid using least squares?
    end_α, end_β, end_Γ = let
        # Precompute lₑ at the starting sigma-points
        lₑ = map(sigma) do
            mahalanobis(sigma[1:L-1], sigma[L], dist)
        end

        # Note, we pass the EuclideanVector thing to set the base type only
        problem = NLLSsolver.NLLSProblem(
            NLLSsolver.EuclideanVector{1, T},
            GVMLeastSquares{L, N, (L * (L-1))÷2, T}
        )
        # end_α
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zero(T)))
        # end_β
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zeros(T, N)...))
        # end_Γ
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zeros(T, S)...))
        # Instantiate the cost-computer
        NLLSsolver.addcost!(problem, GVMLeastSquares(lₑ, endpoints, dist.κ, end_μ, end_A))

        result = NLLSsolver.optimize!(problem, NLLSsolver.NLLSOptions())

        end_α = problem.variables[1][1]
        end_β = problem.variables[2]
        end_Γ_vec = problem.variables[3]

        end_Γ = let
            tmp = zeros(n_Γ, n_Γ)
            tmp[axes(tmp, 1) .>= axes(tmp, 2)'] = end_Γ_vec
            Symmetric(tmp)
        end

        return end_α, end_β, end_Γ

    end

    

end


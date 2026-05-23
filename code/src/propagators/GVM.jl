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
struct GVMSigmaVectors{L, N, T <: Real}
    χ::SMatrix{L, N, T}  # Sigma points matrix
    W::SVector{N, T}     # Weights for each sigma point
end

"""
   b12(κ) = (1 - I₁ / I₀, 1 - I₂ / I₀)

Evaluates modified Bessel functions of the first kind at κ.
Returns the pair efficiently by computing them together.
"""
function b12(κ::T) where {T}
    bessels = Bessels.besseli(0:2, κ)
    return (T(1) - bessels[2] / bessels[1], T(1) - bessels[3] / bessels[1])
end

function GVMSigmaVectors(dist::GaussVonMises{T}) where {T}
    n = length(dist.μ)
    L = n + 1
    N = 2L + 1

    # Canonical sigma vectors depend solely on κ
    b1, b2 = b12(dist.κ)
    ξc = sqrt(T(3))
    ηc = acos(b2 / (T(2) * b1) - T(1))

    wcξ = T(1) / T(6)
    wcη = (b1 * b1) / (T(4) * b1 - b2)
    wc0 = T(1) - T(2) * wcη - T(2) * n * wcξ

    # Generate canonical sigma vectors
    χ = begin
        χbuf = MMatrix{L, N, T}(undef)
        χbuf[:, 1] = SA[zero(dist.μ)..., T(0)]   # χ00
        χbuf[:, 2] = SA[zero(dist.μ)...,  T(ηc)] # χη0
        χbuf[:, 3] = SA[zero(dist.μ)..., -T(ηc)] # χη1

        # Symmetric vectors around the origin for the Euclidean part
        for i in 1:n
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

Residual struct for NLLSsolver to match the likelihood profiles of the
transformed sigma points against the true target GVM distribution.

Type parameters:
    - S: Number of free parameters in Γ, S = ½L(L-1)
    - L: Problem dimensions, L = n + 1
    - N: Number of sigma points, N = 2L + 1
    - M: Number of euclidean dimensions (n)
"""
struct GVMLeastSquares{S, L, N, M, T} <: NLLSsolver.AbstractResidual
    lₑ::SVector{N, T}
    end_σ::SMatrix{L, N, T}
    κ::T
    end_μ::SVector{M, T}
    end_A::LowerTriangular{T, SMatrix{M, M, T}}
end

# NLLSsolver Boilerplate
Base.eltype(::GVMLeastSquares{S, L, N, M, T}) where {S, L, N, M, T} = T
NLLSsolver.ndeps(::GVMLeastSquares) = static(3)
NLLSsolver.nres(::GVMLeastSquares) = static(1)
NLLSsolver.varindices(::GVMLeastSquares) = SVector(1, 2, 3)

function NLLSsolver.getvars(::GVMLeastSquares{S, L, N, M, T}, vars::Vector) where {S, L, N, M, T}
    return (
        vars[1]::NLLSsolver.EuclideanVector{1, T},  # end_α
        vars[2]::NLLSsolver.EuclideanVector{M, T},  # end_β
        vars[3]::NLLSsolver.EuclideanVector{S, T},  # end_Γ (flattened)
    )
end

function NLLSsolver.computeresidual(res::GVMLeastSquares{S, L, N, M, T}, end_α, end_β, end_Γ_vec) where {S, L, N, M, T}
    # 1. Dynamically find the incoming type (Float64 during evaluation, Dual during AD)
    T_AD = promote_type(eltype(end_α), eltype(end_β), eltype(end_Γ_vec))
    n_Γ = L - 1
    
    end_Γ = begin
        tmp = zeros(T_AD, n_Γ, n_Γ)
        tmp[triu(ones(Bool, n_Γ, n_Γ))] = end_Γ_vec
        Symmetric(tmp)
    end

    # 2. Lift the background constants to match the incoming AD type
    μ_ad = T_AD.(res.end_μ)
    A_ad = LowerTriangular(T_AD.(res.end_A))
    κ_ad = T_AD(res.κ)

    # 3. Ensure end_α is treated as a scalar regardless of how NLLSsolver packs it
    α_scalar = end_α isa AbstractArray ? end_α[1] : end_α

    end_dist = GaussVonMises(μ_ad, α_scalar, end_β, end_Γ, κ_ad, A=A_ad)
    
    return sum(zip(res.lₑ, eachcol(res.end_σ))) do (lₑᵢ, end_σᵢ)
        lₐᵢ = mahalanobis(end_σᵢ[1:n_Γ], end_σᵢ[L], end_dist)
        r = T_AD(lₑᵢ) - lₐᵢ
        return r * r
    end
end

function gvm_propagate(f, dist::GaussVonMises{T}) where {T}
    n = length(dist.μ)
    L = n + 1
    N = 2L + 1
    S = (L * (L - 1)) ÷ 2

    # Step 1: Generate and propagate sigma vectors
    sigma = GVMSigmaVectors(dist)
    endpoints = Vector{SVector{L, T}}(undef, N)
    
    Threads.@threads for i in 1:N
        endpoints[i] = f(sigma.χ[:, i])
    end

    # Step 2: GVM quadrature for the Euclidean part
    end_μ, end_P, end_A = let
        euclid_index = [1; 4:N]
        out_μ = sum(endpoints .* sigma.W)[1:n]
        
        dx = reduce(hcat, endpoints)[1:n, euclid_index] .- out_μ
        euclid_W = sigma.W[euclid_index]
        
        out_P = nearest_pd_matrix(dx * Diagonal(euclid_W) * dx')
        out_A = cholesky(Symmetric(out_P)).L
        out_μ, out_P, out_A
    end

    # Step 3: Initial estimate of α, β, and Γ using ForwardDiff
    end_α, end_β, end_Γ = let
        central_point = sigma.χ[:, 1]
        
        δf = ForwardDiff.jacobian(f, central_point)
        δₓfx = δf[1:n, 1:n]
        δₓfα = δf[L, 1:n]

        δ²fα = ForwardDiff.hessian(x -> f(x)[L], central_point)
        δ²ₓfα = δ²fα[1:n, 1:n]

        canon_δₓfx = inv(end_A) * δₓfx * dist.A
        canon_δ²ₓfα = dist.A' * δ²ₓfα * dist.A 

        Δβ = dist.A' * δₓfα
        β_orig = dist.β + Δβ

        α_est = endpoints[1][L]
        β_est = inv(end_A) * δₓfx * dist.A * β_orig
        Γ_est = canon_δₓfx * (dist.Γ + canon_δ²ₓfα) * canon_δₓfx'

        α_est, β_est, Γ_est
    end

    # Step 4: Run Least Squares Refinement in Log Space
    end_α, end_β, end_Γ = let
        lₑ = [mahalanobis(s[1:n], s[L], dist) for s in eachcol(sigma.χ)]

        # Bug fixed: Adjusted type parameter layouts to explicitly prevent inversion crashes
        problem = NLLSsolver.NLLSProblem(Any, GVMLeastSquares{S, L, N, n, T})
        
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zero(T)))
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zeros(T, n)...))
        NLLSsolver.addvariable!(problem, NLLSsolver.EuclideanVector(zeros(T, S)...))
        
        cost_function = GVMLeastSquares{S, L, N, n, T}(
            SVector(lₑ), 
            SMatrix{L, N, T}(reduce(hcat, endpoints)), 
            dist.κ, 
            SVector{n}(end_μ), 
            LowerTriangular(SMatrix{n, n}(end_A))
        )
        NLLSsolver.addcost!(problem, cost_function)

        NLLSsolver.optimize!(problem, NLLSsolver.NLLSOptions())

        α_final = problem.variables[1][1]
        β_final = problem.variables[2]
        Γ_vec   = problem.variables[3]

        Γ_final = begin
            tmp = zeros(T, n, n)
            tmp[triu(ones(Bool, n, n))] = Γ_vec
            Symmetric(tmp)
        end

        α_final, β_final, Γ_final
    end

    return GaussVonMises(end_μ, end_α, end_β, end_Γ, dist.κ, A=end_A)
end

"""
    run_gvm(p::ForceModel, dist::GaussVonMises, Δt; reltol, abstol)

Runs the GVM propagation for a given orbital distribution over time-interval Δt.
"""
function run_gvm(p::ForceModel, dist::GaussVonMises{T, V}, Δt, reltol=1e-10, abstol=1e-10) where {T, V}
    return gvm_propagate(dist) do v
        propagate_orbit(p, v, Δt, reltol=reltol, abstol=abstol)
    end
end

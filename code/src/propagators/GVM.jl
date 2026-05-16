# GVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Gauss Von-Mises propagator, essentially the unscented transform adapted to GVM distribution.
# Implemented according to 
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


function gvm_propagate(
    f,
    dist::GaussVonMises{T,V,M}
) where {T,V,M}

    n = length(dist.μ)
    L = n + 1
    N = 2 * L + 1

    # Step 1: generate the sigma vectors and propagate them
    sigma = GVMSigmaVectors(dist)
    endpoints = Vector{SVector{L,T}}(undef, N)

    Threads.@threads for i in 1:N
        endpoints[i] = f(sigma.χ[:, i])
    end

    # Step 2: GVM cuadrature on the resulting sigma-vectors for the
    #         euclidean part of the distribution
    end_μ = sum(endpoints .* sigma.W)
    dx = reduce(hcat, endpoints) .- end_μ
    end_P = nearest_pd_matrix(dx * Diagonal(sigma.W) * dx')
    end_A = cholesky(Symmetric(end_P)).L

    # Step 3: Estimate hα, hβ, hΓ using derivatives of f
    # RESEARCH: It would be wise to experiment using the derivatives that
    # can be inferred from the endpoints instead of a full forward diff of f.

    # δf[i, j] = δfᵢ/δxⱼ, so we have...
    δf = ForwardDiff.jacobian(f, sigma.N[:, 1])
    δₓfx = δf[1:(L-1), 1:(L-1)]
    δₓfα = δf[L, 1:(L-1)]

    # δ²fθ[i, j] = δ²f / δxᵢδxⱼ
    δ²fα = ForwardDiff.hessian(x -> f(x)[L], sigma.N[:, 1])
    δ²ₓfα = δ²fα[1:(L-1), 1:(L-1)]

    # end_α = f_α(μ, α), that's one of the points we propagate so trivial
    end_α = endpoints[1][L]

    # Total uncertainty in angular coordinates is contributed by original uncertainty β and
    # the new linear contribution from f, after we bring it into canonical coordinates
    # Note that δfₓfα is a 1-form so it transforms by the transpose!
    Δβ = dist.A' * δₓfα

    # Both dist.β and Δβ are 1-forms that act on the original canonical space, thus they can be
    # added together just fine.
    end_β_original_coords = (dist.β + Δβ)

    # Now, we wish to find the transformation that goes from 1-forms in the canonical space before "f", to
    # 1-forms in the canonical space after "f", as a linear approximation of course.
    # We note that δₓfx * ◌ maps vectors from real space before x to real space after x, thus
    # δₓfx * A * ◌ maps a vector in canonical space before x, to a vector in real space after x,
    # and thus (end_A)⁻¹ * δₓfx * A * ◌ maps a vector in canonical space before x, to a vector in canonical space after x.
    canon_δₓfx = inv(end_A) * δₓfx * dist.A

    end_β = inv(end_A) * δₓfx * dist.A * (dist.β + Δβ)
    # Intuitively:
    #   δₓfx represents the linear approximation of f on the euclidean distribution
    #   end_A⁻¹ finally brings everything to the new canonical coordinates

    canon_δₓfx = inv(end_A) * δₓfx * dist.A
    # Intuitively:
    #    let m1 = δₓfx * dist.A * ◌ : δx (canonical coordinates) -> δ̂x (real coordinates) 
    #    let m2 = inv(end_A) * ◌ : δ̂x (real coordinates) -> δ̂x (canonical coordinates) by inverse definition of A
    # then
    #    m * ◌ = m2 ∘ m1 : δx (canonical coordinates ) -> δ̂x (canonical coordinates)
    # i.e. m: ℝⁿ -> ℝⁿ maps a perturbation δx, written in canonical coordinates, to the perturbation in f(x), written in
    #      canonical coordinates of the resulting Gauissian. It's precisely the canonicalized Jacobian.

    # Now lets consider the form Γ: ℝⁿ -> ℝⁿ -> ℝ, how does it transform under f? Note that Γ acts on canonical vectors.
    # To understand it, consider another bilinear symmetric form, the Hessian, interpreting it as a
    # quadratic form Q(δ²ₓfα, ◌) = ◌ᵀ δ²ₓfα ◌ : ℝⁿ -> ℝ, that maps each perturbation δx(real coordinates) to a increase δα.
    # The mapping is NOT a 1-form as it's not linear, but we can derive an expression that allows it to act on δx (canonical coordinates).
    # To do so, consider Q(δ²ₓfα, A δx) = (A δx)ᵀ (δ²ₓfα) (A δx) = (δxᵀ Aᵀ) (δ²ₓfα) (A δx),
    # by associativity of matrix multiplication, we can write δxᵀ (Aᵀ δ²ₓfα A) δx, thus this is precisely the canonicalized Hessian:
    canon_δ²ₓfα = dist.A' + δ²ₓfα  * dist.A 

    # Thus assuming linearized behaviour, ΔΓ = canon_δ²ₓfα, but we must express this
    # "end_Γ" = dist.Γ + ΔΓ in the new coordinates, which is achieved by action of the canonicalized Jacobian 
    end_Γ = canon_δₓfx * (dist.Γ + canon_δ²ₓfα) * canon_δₓfx'


    # Step 4: Solve least square problem to refine hα, hβ, hΓ
    # RESEARCH: If we set Γ=0, do we drop the need for least squares?
    # TODO

end


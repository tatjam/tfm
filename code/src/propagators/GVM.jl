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


function gvm_propagate(f, dist::GaussVonMises{T,V,M}) where {T,V,M}
    n = length(dist.μ)
    L = n + 1
    N = 2 * L + 1

    # Step 1: generate the sigma vectors and propagate them
    sigma = GVMSigmaVectors(dist)
    endpoints = Vector{SVector{L, T}}(undef, N)

    Threads.@threads for i in 1:N
        endpoints[i] = f(sigma.χ[:, i])
    end

    # Step 2: GVM cuadrature on the resulting sigma-vectors

    # Step 3: Estimate hα, hβ, hΓ using derivatives of f
    # RESEARCH: It would be wise to experiment using the derivatives that
    # can be inferred from the endpoints instead of a full forward diff of f.

    # Step 4: Solve least square problem to refine hα, hβ, hΓ
    # RESEARCH: If we set Γ=0, do we drop the need for least squares?

end


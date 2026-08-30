# PolyChaos.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Polynomial chaos expansion based propagator, using the non-intrusive approach.
# Coordinates may either be euclidean or angular. The main idea is our
# chaos inputs are independent, thus the polynomial chaos expansion
# is just the product of the polynomials of both (truncated to maximum
# order). 


"""
A polynomial basis orthogonal to some distribution
"""
abstract type AbstractPCEBasis end

"""
   nvars(b::AbstractPCEBasis)

Return the number of variables in the basis 'b', i.e. the dimensionality of ξ 
"""
function nvars end

"""
    multi_index(b::AbstractPCEBasis)

Returns a vector of tuples, one per basis function, in the same order as used by
eval_basis. Each tuple contains as many entries as variables in the basis, each entry
thus indicating the polynomial degree used for said variable in said basis function.

"""
function multi_index end

"""
Number of basis functions in the basis. 
"""
Base.length(b::AbstractPCEBasis) = length(multi_index(b))

"""
Resulting type of evaluating the basis functions. 
"""
Base.eltype(b::AbstractPCEBasis) = error("Not implemented for type $(typeof(b))")

"""
    eval_basis!(out::AbstractVector, b::AbstractPCEBasis, ξ)

Evaluates all of the basis functions at a given point, setting
    out[i] = Ψ_{α_i}(ξ)
for each α in the multi index, and returning a reference to out.

Note the indices represent order + 1, i.e. out[1] is the order 0, and is
always equal to 1.0 as we use normalized basis.
"""
function eval_basis! end

eval_basis(b::AbstractPCEBasis, ξ) = eval_basis!(Vector{eltype(b)}(undef, length(b)), b, ξ)

abstract type AbstractPCEQuadrature end

"""
    nodes(q::AbstractPCEQuadrature)

Returns a vector of quadrature nodes (values of ξ) 
"""
function nodes end

"""
    weights(q::AbstractPCEQuadrature)

Returns a weight for each of the quadrature nodes
"""
function weights end

"""
   quadrature(b::AbstractPCEBasis) 

Returns a quadrature consistent with the orthogonality measure of the basis.
"""
function quadrature end

"""
   galerkin(f, b::AbstractPCEBasis, q::AbstractPCEQuadrature = quadrature(b))

Obtains the coefficients for the basis b by Galerkin projection of f onto it. The Galerkin
projection means that we choose as a representation of f, such that the difference of the
representation of f to f itself is orthogonal to our basis. We say
projection because geometrically we are finding the projection of f into the subspace defined
by our basis. Refer to fᵣ as

    fᵣ(ξ) = ∑ᵦ cᵦ Ψᵦ(ξ)

then, Galerkin imposes that f - fᵣ is orthogonal to our basis, which means it's orthogonal to every
individual basis function γ,

    <f - ∑ᵦ cᵦ Ψᵦ, Ψᵧ> = 0

by linearity of the expected value, we get

    ∑ᵦ cᵦ <Ψᵦ, Ψᵧ> = <f, Ψᵧ>

This can be made efficient via using a orthonormal base, as then <Ψᵦ, Ψᵧ> = 1 iff β = γ, 0 otherwise,
and thus we get a series of equations for each β

    cᵦ = <f, Ψᵦ>

Finally, we evaluate the expected value via the quadrature rule, which is chosen so it's highly accurate for
our distribution (in fact, ideally should be exact if f is a linear combination of the basis functions). Note that
if high performance is desired, precomputation of the quadrature will save some effort.
"""
function galerkin(f, b::AbstractPCEBasis, q::AbstractPCEQuadrature=quadrature(b))
    ns, ws = nodes(q), weights(q)
    T = eltype(b)
    n = length(b)

    probe_idx = 1
    probe = f(ns[probe_idx])
    m = length(probe)

    # TODO: This could waste a whole lot of memory if maxthreadis is very big!
    partials = [zeros(T, n, m) for _ in 1:Threads.maxthreadid()]
    bufs = [Vector{T}(undef, n) for _ in 1:Threads.maxthreadid()]

    Threads.@threads :static for k in eachindex(ns)
        tid = Threads.threadid()
        buf = bufs[tid]
        eval_basis!(buf, b, ns[k])
        if k == probe_idx
            fk = probe
        else
            fk = f(ns[k])
        end
        partials[tid] .+= ws[k] .* (conj.(buf) * transpose(fk isa Number ? [fk] : fk))
    end

    return sum(partials)
end

"""
A **normalized** polynomial basis orthogonal to some distribution on the reals,
OPRL meaning Orthogonal Polynomial on the Real Line. 

The recurrence relation is compactly stored,
jacobi[i] yields the Favard terms cᵢ₋₁ and √dᵢ in that order, for the relation
  yₙ₊₁= (x - cₙ) yₙ - dₙ yₙ₋₁

which we normalize to 
  pₙ₊₁ = yₙ / sqrt(hₙ) = ((x - cₙ) pₙ - sqrt(dₙ) pₙ₋₁) / sqrt(dₙ₊₁)

Note that, to compute pₙ, {cₙ₋₁, dₙ₋₁ and dₙ} are needed. This justifies storing
the coefficients as follows:

  jacobi[1] = (c₀, √d₁)
  jacobi[2] = (c₁, √d₂)
  ...
  jacobi[n=end] = (cₙ₋₁, √dₙ)

Finally, note that the order N polynomial is NOT usable for PCE, so evaluating the basis
and similar gives results only up to order N-1. This node is instead used only for the quadrature.
"""
struct OprlBasis{T<:AbstractFloat,N} <: AbstractPCEBasis
    jacobi::SVector{N,Tuple{T,T}}
end

# N-1 is the maximum order usable in the basis, plus the 0 order, gives N
Base.length(b::OprlBasis{T, N}) where {T, N} = N 
Base.eltype(b::OprlBasis{T, N}) where {T, N} = T

function OprlBasis(terms::AbstractVector{Tuple{T,T}}) where {T<:AbstractFloat}
    N = length(terms)
    return OprlBasis{T,N}(SVector{N,Tuple{T,T}}(terms))
end

nvars(b::OprlBasis) = 1
multi_index(_b::OprlBasis{T,N}) where {T,N} = 0:N

function eval_basis!(out::AbstractVector{T}, b::OprlBasis{T,N}, ξ) where {T,N}
    # The Favard form is just
    #   yₙ₊₁= (x - cₙ) yₙ - dₙ yₙ₋₁
    # Let hₙ = <yₙ, yₙ>, we wish to construct instead the sequence
    #   pₙ = yₙ / sqrt(hₙ)
    # which is orthonormal. Now, note that Favard sequences satisfy
    #   <yₙ, yₙ> = dₙ <yₙ₋₁, yₙ₋₁>
    # thus hₙ = dₙ hₙ₋₁
    # Now consider pₙ₊₁
    #  pₙ₊₁= yₙ₊₁/ (sqrt(dₙ₊₁hₙ)
    #          = [(x - cₙ) yₙ / sqrt(hₙ) - dₙ yₙ₋₁ / sqrt(hₙ)] / sqrt(dₙ₊₁)
    #          = [(x - cₙ) pₙ            - dₙ yₙ₋₁ / sqrt(dₙ hₙ₋₁)] / sqrt(dₙ₊₁)
    #          = [(x - cₙ) pₙ - sqrt(dₙ) pₙ₋₁] / sqrt(dₙ₊₁)

    # y₀ = p₀ = 1
    out[1] = one(T)

    # y₁ = (x - c₀) y₀, thus
    # sqrt(d₁) p₁ = (x - c₀) p₀
    if N >= 1
        c0, sqd1 = b.jacobi[1]
        out[2] = ((ξ - c0) * out[1]) / sqd1
    end

    # yₙ₊₁= (x - cₙ) yₙ - dₙ yₙ₋₁, thus
    # sqrt(dₙ₊₁ pₙ₊₁= (x - cₙ) pₙ - sqrt(dₙ) pₙ₋₁
    # Note the last element is ignored, that's only for quadrature!
    for n in 1:(N-2)
        cn, sqdnp1 = b.jacobi[n+1]
        _, sqdn = b.jacobi[n]

        # the +1 due to 1 based indexing
        out[(n+1)+1] = ((ξ - cn) * out[n+1] - sqdn * out[(n-1)+1]) / sqdnp1
    end

    return out
end

"""
Basis orthogonal to exp(-0.5x²) in the real line (Probabilist's Hermite polynomials) for
orders up to N (included) usable for PCE (thus the basis stores one more)
"""
function hermite_basis(::Type{T}, N::Int) where {T<:AbstractFloat}
    # The Hermite generator is
    #   yₙ₊₁ = x yₙ - n yₙ₋₁
    # thus cₙ = 0, dₙ = n
    jacobi = [(zero(T), sqrt(T(n))) for n in 1:N+1]
    # @Main.infiltrate

    return OprlBasis(jacobi)
end


struct OprlQuadrature{T<:AbstractFloat,N} <: AbstractPCEQuadrature 
    nodes::SVector{N,T}
    weights::SVector{N,T}
end

function quadrature(b::OprlBasis{T,N}) where {T,N}
    # Note the Jacobi matrix is written by Golub and Welsch as
    #   yₙ = (aₙx + bₙ)yₙ₋₁ - cₙyₙ₋₂
    #   [-b₁/a₁   1/a₁                   ...
    #   [c₂/a₂    -b₂/a₂   1/a₂          ...
    #   [         c₃/a₃    -b₃/a₃   1/a₃ ...
    #
    # But our recurence relation is
    #  pₙ₊₁= [(x - cₙ) pₙ - sqrt(dₙ) pₙ₋₁] / sqrt(dₙ₊₁)
    #
    # Identifying terms (see writeup), we write the Jacobi matrix as
    #   [c₀   √d₁       ...
    #   [√d₁  c₁   √d₂  ...
    #   [     √d₂  c₃   ...

    # The diagonal is all the first terms of jacobi, the off-diagonal the second ones,
    # note that as expected there's one "unused" element in the jacobi array
    jacobimat = SymTridiagonal(Vector(first.(b.jacobi)), Vector(last.(b.jacobi[SOneTo(N - 1)])))

    eigdecomp = eigen(jacobimat)
    eigval = eigdecomp.values
    eigvec = eigdecomp.vectors

    # By Golub and Welsch (Calculation of Gauss quadrature rules, 1969), the
    # nodes of the quadrature are the eigenvalues of the Jacobi matrix:
    nodes = eigval
    # and the weights are simply the first coordinates of each eigen vector squared
    # (multiplied by the total mass of the distribution, but by definition that's 1)
    weights = eigvec[1, :] .^ 2

    return OprlQuadrature(SVector{N}(nodes), SVector{N}(weights))
end


nodes(q::OprlQuadrature{T, N}) where {T, N} = q.nodes
weights(q::OprlQuadrature{T, N}) where {T, N} = q.weights


"""
A polynomial basis orthogonal to some distribution on the unit circle,
OPUC meaning Orthogonal Polynomial on the Unit Circle
"""
struct OpucBasis <: AbstractPCEBasis

end


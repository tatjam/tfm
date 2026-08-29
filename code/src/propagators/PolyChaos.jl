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
our distribution (in fact, ideally should be exact if f is a linear combination of the basis functions).

NOTE: Assumes f has as many outputs as stochastic inputs, could be easily adapted to be more general
"""
function galerkin(f, b::AbstractPCEBasis, q::AbstractPCEQuadrature = quadrature(b))
    ns, ws = nodes(q), weights(q)
    T = eltype(b)
    n = length(b)
    m = nvars(b)

    partials = [zeros(T, n, m) for _ in 1:nthreads()]
    bufs = [Vector{T}(undef, n) for _ in 1:nthreads()]


    # Quadrature rule just means we replace integration by a sum over
    # nodes with weights, so loop over each node. The :static prevents
    # threadid from changing unexpectedly!
    @Threads.threads :static for k in eachindex(ns)
        tid = Threads.threadid()
        buf = bufs[tid]
        eval_basis!(buf, b, ns[k])

        partials[tid] .+= ws[k] .* (conj.(buf) * transpose(f(ns[k])))
    end

    return sum(partials)
end

"""
A polynomial basis orthogonal to some distribution on the reals,
OPRL meaning Orthogonal Polynomial on the Real Line. 

We store the Favard's recurrence terms, as those represent enough
polynomials basis for our purposes.
"""
struct OprlBasis{T<:AbstractFloat} <: Basis
    degree::Int
    favard::{Tuple{T, T}}
end

nvars(b::OprlBasis) = 1
multi_index(b::OprlBasis) = 0:b.degree
function eval_basis(b::OprlBasis{T}, ξ::T) where T 
    
end


"""
A polynomial basis orthogonal to some distribution on the unit circle,
OPUC meaning Orthogonal Polynomial on the Unit Circle
"""
struct OpucBasis <: Basis

end

function run_poly_chaos(p::ForceModel, )


end

# PolyChaos.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Polynomial chaos expansion based propagator, using the non-intrusive approach.
# Coordinates may either be euclidean or angular. The main idea is our
using Distributions: degrees
# chaos inputs are independent, thus the polynomial chaos expansion
# is just the product of the polynomials of both (truncated to maximum
# order). 


"""
A polynomial basis orthogonal to some distribution, if ξ is our vector
of random variables, this represents a sum of a number of polynomial functions
of ξ, each indexed by a index αᵢ
"""
abstract type AbstractPCEBasis end

"""
   nvars(b::AbstractPCEBasis)

Return the number of variables in the basis 'b', i.e. the dimensionality of ξ 
"""
function nvars end

"""
    multi_index(b::AbstractPCEBasis)

Return a multi_index object, each of its entries is one of the terms in the basis,
so iterating over this yields every possible α.

Note that the returned object could be an "opaque" pointer if the basis
is not a tensor product, this is not neccesarily a conventional multi index.

"""
function multi_index end

"""
    eval_basis(b::Basis, ξ, α)

Evaluates one of the basis functions at a given point.
"""
function eval_basis end

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

function eval_basis(b::OprlBasis{T}, ξ::T, α::Integer) where T 
    
end


"""
A polynomial basis orthogonal to some distribution on the unit circle,
OPUC meaning Orthogonal Polynomial on the Unit Circle
"""
struct OpucBasis <: Basis

end

function run_poly_chaos(p::ForceModel, )


end

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


"""
    gvm(μ, P, α, β, Γ, κ)

Create a Gauss Von Mises distribution for n Euclidean dimensions and a single angular dimension. The distribution is defined as N(μ, P) × VM(α + β'z + 1/2 z'Γz, κ). z = A^(-1)(x - μ). Where A is the lower triangular Cholesky decomposition of P, and:

     μ is a n dimensional vector representing the mean of the Euclidean elements
     P is a n×n matrix representing the covariance of the Euclidean elements. It's thus symmetric and positive definite.
     α is a real number representing the "angular mean" of the angular element
     β is a n dimensional vector capturing the linear correlation between the euclidean elements and the angular element
     Γ is a n×n matrix used as a quadratic form to tune the shape of the distribution.
     κ is a real number representing the angular concentration of the distribution.

The distribution assigns probabilities to a random tuple (x, θ), where x is the euclidean random variable and θ is the angular random variable.
    
"""
function gvm(μ, P, α, β, Γ, κ)

end

# STM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# State Transition Matrix Gaussian uncertainty propagation.
# We use automatic differentation to build the state transition matrix, and use
# this matrix to update the statistical distribution

using ForwardDiff
using DiffResults

"""
    run_stm(p::ForceModel, dist::Distribution, Δt; reltol=1e-10, abstol=1e-10) 

Propagates the dynamics of the force model and builds the state transition matrix (STM) by
automatic differentiation. This STM, alongside the propagated mean,
is used to update the distribution, which is returned.
"""
function run_stm(
    p::ForceModel,
    μ::AbstractVector{T},
    P::AbstractMatrix{T},
    Δt;
    reltol=1e-10,
    abstol=1e-10,
) where {T}

    dr = DiffResults.JacobianResult(μ)
    dr = ForwardDiff.jacobian!(dr, x -> propagate_orbit(p, x, Δt, reltol=reltol, abstol=abstol), μ)

    μ_prop = DiffResults.value(dr)
    # Note, this is not really the jacobian as the entire integration is performed, this
    # is directly the state transition matrix
    ϕ = DiffResults.jacobian(dr)

    Pend = nearest_pd_matrix(ϕ * P * ϕ')
    return MvNormal(μ_prop, Pend)
end

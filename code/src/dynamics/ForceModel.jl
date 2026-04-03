# ForceModel.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Generalized forces to make the propagator easily configurable,
# using a compile-time Tuple for high performance dispatch.
#
# Both "Newton" propagation (Chapter 8 of Vallado) and "Gauss" (Chapter 9 of Vallado) 
# propagation are possible with the same data structures.

"""
    TwoBodyForce(μ)

Force due to a central body at the origin of the coordinate system, with gravitational
parameter `μ`.
"""
struct TwoBodyForce
    μ::Float64
end

"""
    acceleration(f::TwoBodyForce, r, _, _)

Newton's law of universal gravitation
"""
function acceleration(f::TwoBodyForce, r, _v, _t)
    -f.μ / norm(r)^3 * r
end

"""
    param_variation(fm::TwoBodyForce, p, f, g, h, k, L, t)

    "A set of modified equinoctial orbit elements", Walker et al 1985, formula 9
    with all perturbations set to 0

"""
function param_variation(fm::TwoBodyForce, p, f, g, _h, _k, L, _t)
    w = 1 + f * cos(L) + g * sin(L)

    # Central force merely makes L increase
    dL = sqrt(fm.μ * p) * (w / p)^2

    return SA[0, 0, 0, 0, 0, dL]
end

"""
    J2Force(μ, R, J2)

J2 perturbation due to a central body at the origin of the coordinate system
"""
struct J2Force
    μ::Float64
    R::Float64
    J2::Float64
end

"""
    acceleration(f::J2Force, r, v, t)

J2 newtonian perturbation, Vallado page 594 formula.
"""
function acceleration(f::J2Force, r, _v, _t)
    common = -3.0 * f.J2 * f.μ * f.R^2 / (2.0 * norm(r)^5)
    zrel = 5 * r[3]^2 / norm(r)^2
    xy_term = 1.0 - zrel
    z_term = 3.0 - zrel
    return common * (r .* SA[xy_term, xy_term, z_term])
end

"""
    param_variation(fm::J2Force, p, f, g, h, k, L, t)

    "A set of modified equinoctial orbit elements", Walker et al 1985, formula 8
    and formula 11 particularized for J2 only

"""
function param_variation(fm::J2Force, p, f, g, h, k, L, _t)
    sinL, cosL = sincos(L)
    w = 1 + f * cosL + g * sinL
    s = sqrt(1 + h^2 + k^2)
    r = p / w
    sinϕ = 2 * (h * sinL - k * cosL) / s^2
    # P_2(sinϕ)
    Pn = 0.5 * (3 * sinϕ^2 - 1)
    # P_2'(sinϕ)
    dPn = 3 * sinϕ
    cterm = fm.J2 * (fm.R / r)^2

    # Potential gradient (formula 11, with TYPO fix!)
    dRdp = 3 * fm.μ / (w * r^2) * cterm * Pn
    # TYPO on formula 11, dRdf is missing, it's the symmetrical to dRdg
    # See ERRATA for the Walker et al paper (https://link.springer.com/article/10.1007/BF01238929)
    dRdf = -3 * fm.μ * cosL / (w * r) * cterm * Pn

    dRdg = -3 * fm.μ * sinL / (w * r) * cterm * Pn
    dRdh = -2 * fm.μ / (r * s^4) * ((1 - h^2 + k^2) * sinL + 2 * h * k * cosL) * cterm * dPn
    dRdk = 2 * fm.μ / (r * s^4) * ((1 + h^2 - k^2) * cosL + 2 * h * k * sinL) * cterm * dPn
    dRdL1 = -2 * fm.μ / (r * s^2) * (h * cosL + k * sinL) * cterm * dPn
    dRdL2 = -3 * fm.μ / (r * w) * (g * cosL - f * sinL) * cterm * Pn
    dRdL = dRdL1 + dRdL2

    # Gauss equations (formula 8)
    sqrtμp = sqrt(fm.μ * p)
    dp = 2 * sqrt(p / fm.μ) * (-g * dRdf + f * dRdg + dRdL)
    df1 = 2 * p * g * dRdp - (1 - f^2 - g^2) * dRdg - 0.5 * g * s^2 * (h * dRdh + k * dRdk)
    df2 = (f + (1 + w) * cosL) * dRdL
    df = 1 / sqrtμp * (df1 + df2)
    dg1 = -2 * p * f * dRdp + (1 - f^2 - g^2) * dRdf + 0.5 * f * s^2 * (h * dRdh + k * dRdk)
    dg2 = (g + (1 + w) * sinL) * dRdL
    dg = 1 / sqrtμp * (dg1 + dg2)
    dh = 0.5 * s^2 / sqrtμp * (h * (g * dRdf - f * dRdg - dRdL) - 0.5 * s^2 * dRdk)
    dk = 0.5 * s^2 / sqrtμp * (k * (g * dRdf - f * dRdg - dRdL) + 0.5 * s^2 * dRdh)

    # note the two-body term was removed from dL
    dL = 0.5 * s^2 / sqrtμp * (h * dRdh + k * dRdk)

    return SA[dp, df, dg, dh, dk, dL]
end

"""

    ForceModel{IsNewton}(forces)

A generic force model, that applies forces sequentially. The `forces` argument is a tuple
of various forces, fixed at compile-time for optimal execution.

If IsNewton is true, the ForceModel works with "cartesian" elements:
    [x0, x1, x2, v0, v1, v2]

If IsNewton is false, the ForceModel works with "keplerian" elements. To avoid singularity,
we use the modified equinoctial elements (A Set of Modified Equinoctial Orbit Elements, Walker et al, 1985)
    [p, f, g, h, k, L] 

defined as follows
    p = a (1 - e²)
    f = e cos(ω + Ω)
    g = e sin(ω + Ω)
    h = tan(½i) cos(Ω)
    k = tan(½i) sin(Ω)
    L = Ω + ω + ν
"""
struct ForceModel{F<:Tuple,IsNewton}
    forces::F
end

ForceModel(forces::F, ::Val{IsNewton}) where {F<:Tuple,IsNewton} =
    ForceModel{F,IsNewton}(forces)


"""

    acceleration(fm::ForceModel, r, v, t)

Computes the acceleration due to all forces in the model sequentially
"""
function acceleration(fm::ForceModel{F,true}, r, v, t) where {F}
    a = SA[0.0, 0.0, 0.0]
    for force in fm.forces
        a += acceleration(force, r, v, t)
    end

    return a
end

"""

    param_variation(fm::ForceModel, r, v, t)

Computes the variation in MEE parameters given the force model.
"""
function param_variation(fm::ForceModel{F,false}, u, t) where {F}
    du = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for force in fm.forces
        du += param_variation(force, u..., t)
    end
    return du
end

"""

    force_model(u, p::ForceModel{True}, t)

Propagation equations for a Newtonian ForceModel. Computes du, given u,
the force model, and the time.
"""
function force_model(u, p::ForceModel{F,true}, t) where {F}
    r = SA[u[1], u[2], u[3]]
    v = SA[u[4], u[5], u[6]]
    a = acceleration(p, r, v, t)

    return SA[v[1], v[2], v[3], a[1], a[2], a[3]]
end

"""

    force_model(u, p::ForceModel{False}, t)

Propagation equations for a Keplerian ForceModel. Computes du, given u,
the force model, and the time.
"""
function force_model(u, p::ForceModel{F,false}, t) where {F}
    return param_variation(p, u, t)
end



"""

    propagate_orbit(fm::ForceModel, u0, t; reltol, abstol)

Propagates a ForceModel. Computes du, given u, the force model, and the time.
"""
function propagate_orbit(
    fm::ForceModel,
    u0::SVector{6,<:Real},
    t::Real;
    reltol=1e-10,
    abstol=1e-10
)

    tspan = (zero(t), t)
    prob = ODEProblem(force_model, u0, tspan, fm)
    sol = solve(prob, Tsit5(); reltol, abstol)

    return last(sol.u)
end

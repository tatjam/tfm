# ForceModel.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Generalized forces to make the propagator easily configurable,
# using a compile-time Tuple for high performance dispatch.
#
# Both "Newton" propagation (Chapter 8 of Vallado) and "Gauss" (Chapter 9 of Vallado) 
# propagation are possible with the same data structures.

"""
    J2Force(R, J2)

J2 perturbation due to a central body at the origin of the coordinate system
"""
struct J2Force
    R::Float64
    J2::Float64
end

"""
    acceleration(f::J2Force, r, v, t, μ)

J2 newtonian perturbation, Vallado page 594 formula.
"""
function acceleration(f::J2Force, r, _v, _t, μ)
    # Note this is slightly incorrect if we assume r is on ECI frame, as J2 is defined in
    # earth relative coordinates, and the pole rotates slightly wrt. to our assumed J2000 ECI frame.
    # For the purpose of this work, it doesn't matter, and would slightly impact performance, but,
    # assumign t is relative to J2000 epoch:

    # jd0 = date_to_jd(2000, 1, 1, 1, 1, 1)
    # jd = t / 86400.0 + jd0
    # eci2ecef = r_eci_to_ecef(J2000(), PEF(), jd)
    # r = eci2ecef * r

    common = -3.0 * f.J2 * μ * f.R^2 / (2.0 * norm(r)^5)
    zrel = 5 * r[3]^2 / norm(r)^2
    xy_term = 1.0 - zrel
    z_term = 3.0 - zrel

    # Afterwards, the J2 computation is on correct frame, but needs to be transformed back
    # a = common * (r .* SA[xy_term, xy_term, z_term])
    # a = eci2ecef' * a
    # return a

    return common * (r .* SA[xy_term, xy_term, z_term])
end

"""
    param_variation(fm::J2Force, p, f, g, h, k, L, t, μ)

    "A set of modified equinoctial orbit elements", Walker et al 1985, formula 8
    and formula 11 particularized for J2 only

"""
function param_variation(fm::J2Force, p, f, g, h, k, L, _t, μ)
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
    dRdp = 3 * μ / (w * r^2) * cterm * Pn
    # TYPO on formula 11, dRdf is missing, it's the symmetrical to dRdg
    # See ERRATA for the Walker et al paper (https://link.springer.com/article/10.1007/BF01238929)
    dRdf = -3 * μ * cosL / (w * r) * cterm * Pn

    dRdg = -3 * μ * sinL / (w * r) * cterm * Pn
    dRdh = -2 * μ / (r * s^4) * ((1 - h^2 + k^2) * sinL + 2 * h * k * cosL) * cterm * dPn
    dRdk = 2 * μ / (r * s^4) * ((1 + h^2 - k^2) * cosL + 2 * h * k * sinL) * cterm * dPn
    dRdL1 = -2 * μ / (r * s^2) * (h * cosL + k * sinL) * cterm * dPn
    dRdL2 = -3 * μ / (r * w) * (g * cosL - f * sinL) * cterm * Pn
    dRdL = dRdL1 + dRdL2

    # Gauss equations (formula 8)
    sqrtμp = sqrt(μ * p)
    dp = 2 * sqrt(p / μ) * (-g * dRdf + f * dRdg + dRdL)
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
    EGM96Force(degree, order, jd0)

Gravitational force due to earth under the EGM96 model as exposed by SatelliteToolbox, NOT including the central body
force, only perturbation terms. t0 is the date time for t=0 used for the coordinate transformations to ECEF frame, thus
t is measured in seconds since this epoch.

All coordinates are assumed to be in the J2000 ECI reference frame, to be converted to PEF (Pseudo-Earth fixed, ITRF without polar
motion) for evaluating the gravity model.
Could be trivially adapted to use a more precise earth orientation model, because technically PEF is not correct for EGM96

If degree is negative, the maximum coefficient of EGM96 is used.
If order is negative, a value equal to degree is used.

jd0 is the Julian date for t = 0, for example, use date_to_jd(1986, 6, 19, 21, 35, 0)
"""
struct EGM96Force
    model::AbstractGravityModel{Float64}
    degree::Int32
    order::Int32
    jd0::Float64
    buffer_P::LowerTriangularStorage
    buffer_dP::LowerTriangularStorage

    function EGM96Force(degree, order, jd0)
        model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

        if degree < 0
            degree = GravityModels.maximum_degree(model)
        end
        if order < 0
            order = degree
        end

        buffer_P = LowerTriangularStorage((degree + 1) * (order +1))
        buffer_dP = LowerTriangularStorage((degree + 1) * (order +1))
        new(model, degree, order, jd0, buffer_P, buffer_dP)
    end
end

function egm96_acceleration_eci(f::EGM96Force, r_eci, t)
    # SatelliteToolbox excepts position in ECEF (PEF) frame, but (r, v) are in inertial frame (J2000)
    # Julian date is just time since epoch in days
    jd = t / 86400.0 + f.jd0
    eci2ecef = r_eci_to_ecef(J2000(), PEF(), jd)

    # Note, acc doesn't include rotational terms, and time is expected as seconds, not julian date, since J2000
    acc_ecef = GravityModels.gravitational_acceleration(f.model, eci2ecef * r_eci, jd * 86400.0)
    acc_eci = eci2ecef' * acc_ecef

    # We need to exclude the 2-body term, as it's included by EGM96
    acc_eci += r_eci * GravityModels.gravity_constant(f.model) / norm(r_eci)^3

    return acc_eci
end

"""
    acceleration(f::EGM96Force, r, v, t)

EGM96 newtonian acceleration.
"""
function acceleration(f::EGM96Force, r_eci, _v, t, _μ)
    return egm96_acceleration_eci(f, r_eci, t)
end

"""
    param_variation(fm::EGM96Force, p, f, g, h, k, L, t, μ)

IGRF acceleration, computed in euclidean coordinates and then transformed.
"""
function param_variation(fm::EGM96Force, p, f, g, h, k, L, t, μ)
    euclid_state = mee_to_euclid(p, f, g, h, k, L, μ)
    a_eci = egm96_acceleration_eci(fm, euclid_state[SA[1,2,3]], t)

    csn2eci = get_csn_basis(euclid_state...)
    a_csn = csn2eci' * a_eci

    return csn_acceleration_to_mee(p, f, g, h, k, L, a_csn..., μ)
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
    μ::Float64
    forces::F
end

ForceModel(μ::Float64, forces::F, ::Val{IsNewton}) where {F<:Tuple,IsNewton} =
    ForceModel{F,IsNewton}(μ, forces)


"""

    acceleration(fm::ForceModel, r, v, t)

Computes the acceleration due to all forces in the model sequentially
"""
function acceleration(fm::ForceModel{F,true}, r, v, t) where {F}
    a = SA[0.0, 0.0, 0.0]
    for force in fm.forces
        a += acceleration(force, r, v, t, fm.μ)
    end

    # 2-body term
    a += -fm.μ / norm(r)^3 * r

    return a
end

"""

    param_variation(fm::ForceModel, u, t)

Computes the variation in MEE parameters given the force model.
"""
function param_variation(fm::ForceModel{F,false}, u, t) where {F}
    du = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for force in fm.forces
        du += param_variation(force, u..., t, fm.μ)
    end

    # 2-body term
    w = 1 + u[2] * cos(u[6]) + u[3] * sin(u[6])
    du += SA[0.0, 0.0, 0.0, 0.0, 0.0, sqrt(fm.μ * u[1]) * (w / u[1])^2] 

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
    u0::AbstractVector{<:Real},
    t::Real;
    reltol=1e-10,
    abstol=1e-10
)
    u0_static = SVector{6}(u0)

    tspan = (zero(t), t)
    prob = ODEProblem(force_model, u0_static, tspan, fm)
    sol = solve(prob, Tsit5(); reltol, abstol)

    return last(sol.u)
end

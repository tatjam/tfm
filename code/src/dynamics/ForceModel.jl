# ForceModel.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Generalized forces to make the propagator easily configurable,
# using a compile-time Tuple for high performance dispatch.
#
# Both "Newton" propagation (Chapter 8 of Vallado) and "Gauss" (Chapter 9 of Vallado) 
# propagation are possible with the same data structures.

# Default 0 extra parameters
nextra(::Type{T}) where {T} = 0

# Static indexes for performance
const POS = SA[1, 2, 3]
const VEL = SA[4, 5, 6]

@generated function nextra_total(::Type{F}) where {F<:Tuple}
    n = 0
    for T in F.parameters
        n += nextra(T)
    end
    return n
end

"""
    J2Force(R, J2)

J2 perturbation due to a central body at the origin of the coordinate system
"""
struct J2Force
    R::Float64
    J2::Float64
end

"""
    acceleration(f::J2Force, rv, extra, t, μ)

J2 newtonian perturbation, Vallado page 594 formula.
"""
function acceleration(f::J2Force, rv, _extra, _t, μ)
    r = rv[POS]

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

    acc = common * (r .* SA[xy_term, xy_term, z_term])

    return SA[0.0, 0.0, 0.0, acc...]
end

"""
    param_variation(fm::J2Force, p, f, g, h, k, L, extra, t, μ)

    "A set of modified equinoctial orbit elements", Walker et al 1985, formula 8
    and formula 11 particularized for J2 only

"""
function param_variation(fm::J2Force, p, f, g, h, k, L, _extra, _t, μ)
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

        buffer_P = LowerTriangularStorage((degree + 1) * (order + 1))
        buffer_dP = LowerTriangularStorage((degree + 1) * (order + 1))
        new(model, degree, order, jd0, buffer_P, buffer_dP)
    end
end

function egm96_acceleration_eci(f::EGM96Force, r_eci, t)
    # SatelliteToolbox excepts position in ECEF (PEF) frame, but (r, v) are in inertial frame (J2000)
    # Julian date is just time since epoch in days
    jd = t / 86400.0 + f.jd0
    eci2ecef = r_eci_to_ecef(J2000(), PEF(), jd)

    # Note, acc doesn't include rotational terms, and time is expected as seconds, not julian date, since J2000
    acc_ecef =
        GravityModels.gravitational_acceleration(f.model, eci2ecef * r_eci, jd * 86400.0)
    acc_eci = eci2ecef' * acc_ecef

    # We need to exclude the 2-body term, as it's included by EGM96
    acc_eci += r_eci * GravityModels.gravity_constant(f.model) / norm(r_eci)^3

    return acc_eci
end

"""
    acceleration(f::EGM96Force, rv, extra, t, μ)

EGM96 newtonian acceleration.
"""
function acceleration(f::EGM96Force, rv_eci, _extra, t, _μ)
    r_eci = rv_eci[POS]
    return SA[0, 0, 0, egm96_acceleration_eci(f, r_eci, t)...]
end

"""
    param_variation(fm::EGM96Force, p, f, g, h, k, L, extra, t, μ)

EGM96 acceleration, computed in euclidean coordinates and then transformed.
"""
function param_variation(fm::EGM96Force, p, f, g, h, k, L, _extra, t, μ)
    euclid_state = mee_to_euclid(p, f, g, h, k, L, μ)
    a_eci = egm96_acceleration_eci(fm, euclid_state[SA[1, 2, 3]], t)

    csn2eci = get_csn_basis(euclid_state...)
    a_csn = csn2eci' * a_eci

    return csn_acceleration_to_mee(p, f, g, h, k, L, a_csn..., μ)
end



"""
   SimpleDragForce(ρ₀, r₀, H)

Simple exponential falloff atmosphere model with density
    ρ(r) = ρ₀ exp(-(r - r₀) / H)
where r is the radius to the earth's center. The atmosphere is assumed
to remain stationary relative to the earth fixed frame for simplicity.

Drag is applied according to the ballistic coefficient β, giving
    a = ½ρv² / β
opposite to the direction of motion.
"""
struct SimpleDragForce{T<:Real}
    # Density (kg/m³) at r₀
    ρ₀::T
    # Radius of reference (m)
    r₀::T
    # Scale height (m)
    H::T
end

# The ballistic coefficient of the body β
nextra(::Type{<:SimpleDragForce}) = 1

"""
    acceleration(f::SimpleDragForce, rv, extra, t, μ)

Drag from the simple atmosphere model.
"""
function acceleration(f::SimpleDragForce, rv, extra, t, _μ)
    r = rv[SA[1, 2, 3]]
    v = rv[SA[4, 5, 6]]
    β = extra[1]

    # Angular velocity of the earth = atmosphere
    ωE = SA[0.0, 0.0, 7.2921159e-5]
    v_rel = v - cross(ωE, r)

    ρ = f.ρ₀ * exp(-(norm(r) - f.r₀) / f.H)
    a = -v_rel * 0.5 * ρ * norm(v_rel) / β
    return SA[0.0, 0.0, 0.0, a..., 0.0]
end

"""
    param_variation(fm::SimpleDragForce, p, f, g, h, k, L, extra, t, μ)

Drag from the simple atmosphere mode
"""
function param_variation(fm::SimpleDragForce, p, f, g, h, k, L, extra, _t, μ)
    euclid_state = mee_to_euclid(p, f, g, h, k, L, μ)
    a_eci = acceleration(fm, euclid_state, extra, _t, μ)[SA[4, 5, 6]]
    csn2eci = get_csn_basis(euclid_state...)
    a_csn = csn2eci' * a_eci

    return SA[csn_acceleration_to_mee(p, f, g, h, k, L, a_csn..., μ)..., 0.0]
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

Either way, force models can have extra elements after the main ones to specify
per-propagation values, as determined by the nextra function. The elements
are stacked one after another in the array in the same order as the force models
are in the tuple, and the extra argument receives only the relevant ones.
"""
struct ForceModel{F<:Tuple,IsNewton,N}
    μ::Float64
    forces::F
end

function ForceModel(μ::Float64, forces::F, ::Val{IsNewton}) where {F<:Tuple,IsNewton}
    N = nextra_total(F)
    return ForceModel{F,IsNewton,N}(μ, forces)
end

# Compile-time sized slice of the extras, explicit eltype so n = 0 works
@inline extras(u, o, ::Val{n}) where {n} =
    SVector{n,eltype(u)}(ntuple(i -> u[o+i-1], Val(n)))

# Recursive to ensure it's optimized away into flat code
@inline add_acceleration!(_du, ::Tuple{}, _u, _o, _t, _μ) = nothing
@inline function add_acceleration!(du, forces::Tuple, u, o, t, μ)
    force = first(forces)
    n = nextra(typeof(force))
    a = acceleration(force, u[SOneTo(6)], extras(u, o, Val(n)), t, μ)
    for i = 1:6
        du[i] += a[i]
    end
    for i = 1:n
        du[o+i-1] += a[6+i]
    end
    add_acceleration!(du, Base.tail(forces), u, o + n, t, μ)
end

# Same as before
@inline add_variation!(_du, ::Tuple{}, _u, _o, _t, _μ) = nothing
@inline function add_variation!(du, forces::Tuple, u, o, t, μ)
    force = first(forces)
    n = nextra(typeof(force))
    a = param_variation(force, u[SOneTo(6)]..., extras(u, o, Val(n)), t, μ)
    for i = 1:6
        du[i] += a[i]
    end
    for i = 1:n
        du[o+i-1] += a[6+i]
    end
    add_variation!(du, Base.tail(forces), u, o + n, t, μ)
end

"""
    acceleration(fm::ForceModel, u, t)

Sum of all force contributions (and extra-parameter derivatives), plus the
two-body term.
"""
function acceleration(fm::ForceModel{F,true,N}, u, t) where {F,N}
    T = promote_type(Float64, eltype(u))
    du = zero(MVector{6 + N,T})
    add_acceleration!(du, fm.forces, u, 7, t, fm.μ)

    # 2-body term
    r = u[POS]
    g = -fm.μ / norm(r)^3 * r
    du[4] += g[1]
    du[5] += g[2]
    du[6] += g[3]

    return SVector(du)
end

"""
    param_variation(fm::ForceModel, u, t)

Variation of the MEE parameters (and extras) given the force model.
"""
function param_variation(fm::ForceModel{F,false,N}, u, t) where {F,N}
    T = promote_type(Float64, eltype(u))
    du = zero(MVector{6 + N,T})
    add_variation!(du, fm.forces, u, 7, t, fm.μ)

    # 2-body term only in L
    w = 1 + u[2] * cos(u[6]) + u[3] * sin(u[6])
    du[6] += sqrt(fm.μ * u[1]) * (w / u[1])^2

    return SVector(du)
end

function force_model(u, p::ForceModel{F,true,N}, t) where {F,N}
    a = MVector(acceleration(p, u, t))
    # dx = v
    a[1] = u[4]
    a[2] = u[5]
    a[3] = u[6]
    return SVector(a)
end

function force_model(u, p::ForceModel{F,false,N}, t) where {F,N}
    return param_variation(p, u, t)
end


"""

    propagate_orbit(fm::ForceModel, u0, t; reltol, abstol)

Propagates a ForceModel. Computes du, given u, the force model, and the time.
"""
function propagate_orbit(
    fm::ForceModel{F,I,N},
    u0::AbstractVector{<:Real},
    t::Real;
    reltol = 1e-10,
    abstol = 1e-10,
) where {F,I,N}
    u0_static = SVector{6+N}(u0)

    tspan = (zero(t), t)
    prob = ODEProblem(force_model, u0_static, tspan, fm)
    sol = solve(prob, Tsit5(); reltol, abstol)

    return last(sol.u)
end

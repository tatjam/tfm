# ForceModel.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Generalized forces to make the propagator easily configurable,
# using a compile-time Tuple for high performance dispatch.

"""Force due to a central body at the origin of the coordinate system"""
struct TwoBodyForce
    μ :: Float64

    TwoBodyForce() = new(SatelliteToolboxBase.GM_EARTH)
end

"""Newton's law of universal gravitation"""
function acceleration(f::TwoBodyForce, r, _v, _t)
    -f.μ / norm(r)^3 * r
end

"""J2 perturbation due to a central body at the origin of the coordinate system"""
struct J2Force
    μ :: Float64
    R :: Float64
    J2 :: Float64 
end

"""J2 perturbation, Vallado page 594""" 
function acceleration(f::J2Force, r, _v, _t)
    common = - 3.0 * f.J2 * f.μ * f.R^2 / (2.0 * norm(r)^2) 
    zrel = 5 * r[3] ^ 2 / norm(r)^2
    xy_term = 1.0 - zrel
    z_term = 3.0 - zrel
    return common * (r .* SA[xy_term, xy_term, z_term])
end


"""A full gravity model from SatelliteToolbox, including the two body force"""
struct GravityModel
    # TODO
end

"""An atmospheric model from SatelliteToolbox"""
struct AtmosphericModel
    # TODO
end

"""A generic force model, that applies forces sequentially"""
struct ForceModel{F <: Tuple}
    forces::F
end

"""Computes the acceleration due to all forces in the model sequentially"""
function acceleration(fm::ForceModel, r, v, t)
    a = SA[0.0, 0.0, 0.0]
    for force in fm.forces
        a += acceleration(force, r, v, t)
    end

    return a
end

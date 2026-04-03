# utils.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Utils for dynamics code 

"""
    kepler_to_mee(a, e, i, ω, Ω, ν)

    Converts from Kepler elements to Modified Equinoctial Elements 
    according to formula (2) of "A set of modified equinoctial orbit elements".

    Returns an array [p, f, g, h, k, L].
"""
function kepler_to_mee(a, e, i, ω, Ω, ν)
    return SA[
        a*(1-e^2),
        e*cos(ω + Ω),
        e*sin(ω + Ω),
        tan(i / 2)*cos(Ω),
        tan(i / 2)*sin(Ω),
        Ω+ω+ν
    ]
end

"""
    kepler_to_mee(k::KeplerianElements)

    Converts from Kepler elements to Modified Equinoctial Elements 
    according to formula (2) of "A set of modified equinoctial orbit elements".

    Returns an array [p, f, g, h, k, L].
"""
function kepler_to_mee(k::KeplerianElements)
    return kepler_to_mee(k.a, k.e, k.i, k.ω, k.Ω, k.f)
end

"""
    mee_to_kepler_array(p, f, g, h, k, L)

    Converts from Modified Equinoctial Elements to Kepler elements, via the inverses
    derived in the Mathematica annex

    Returns an array [a, e, i, ω, Ω, ν]

"""
function mee_to_kepler(p, f, g, h, k, L)
    e = sqrt(f^2 + g^2)
    Ω = atan(k, h)
    ω = atan(g, f) - Ω
    return SA[
        p/(1-e^2),
        e,
        2*atan(sqrt(h^2 + k^2)),
        ω,
        Ω,
        L-ω-Ω
    ]
end

"""
    kepler_to_array(ke::KeplerianElements)

    Returns the KeplerianElements object as a array [a, e, i, ω, Ω, ν]

"""
function kepler_to_array(ke::KeplerianElements)
    return SA[ke.a, ke.e, ke.i, ke.ω, ke.Ω, ke.f]
end

"""
    kepler_to_euclid(a, e, i, ω, Ω, ν, μ)

    Converts kepler elements to euclidean array [x1, x2, x3, v1, v2, v3]

    Implementation is copied from SatelliteToolbox, only big change is that
    the fixed GM_EARTH constant is now passable externally for coherency. 
"""
function kepler_to_euclid(a, e, i, ω, Ω, ν, μ)
    !(0 <= e < 1) && throw(ArgumentError("Eccentricity must be in the interval [0,1)."))

    sin_ν, cos_ν = sincos(ν)
    e2 = e * e

    r = a * (1 - e2) / (1 + e * cos_ν)

    r_o = SA[r * cos_ν, r * sin_ν, 0]

    n₀  = √(μ / a^3)
    v_o = n₀ * a / sqrt(1 - e2) * SA[-sin_ν, e + cos_ν, 0]

    D_i_o = angle_to_dcm(-ω, -i, -Ω, :ZXZ)

    r_i = D_i_o * r_o
    v_i = D_i_o * v_o

    return SA[r_i..., v_i...]
end

"""
   isapprox_angle(a, b; atol=1e-8)

   Compare two angles knowing they live in the circle
"""
isapprox_angle(a, b; atol=1e-8) = abs(rem2pi(a - b, RoundNearest)) < atol


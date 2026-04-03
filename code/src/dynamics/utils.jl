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

isapprox_angle(a, b; atol=1e-8) = abs(rem2pi(a - b, RoundNearest)) < atol


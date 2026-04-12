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
    euclid_to_kepler(x1, x2, x3, v1, v2, v3, μ)

Converts Euclidean elements to Kepler elements [a, e, i, ω, Ω, ν]

Implementation is copied from SatelliteToolbox, with the fixed GM_EARTH removed
"""
function euclid_to_kepler(x1::T, x2::T, x3::T, v1::T, v2::T, v3::T, μ::T) where {T}
    @inbounds begin
        sr_i = SVector{3, T}(x1, x2, x3)
        sv_i = SVector{3, T}(v1, v2, v3)

        r² = dot(sr_i, sr_i)
        v² = dot(sv_i, sv_i)
        r  = sqrt(r²)
        v  = sqrt(v²)
        rv = dot(sr_i, sv_i)

        h_i = sr_i × sv_i
        h   = norm(h_i)

        # Vector that points to the right ascension of the ascending node (RAAN).
        n_i = SVector{3}(0, 0, 1) × h_i
        n   = norm(n_i)

        e_i = ((v² - μ / r) * sr_i - rv * sv_i ) / μ

        # Orbit energy.
        ξ = v² / 2 - μ / r

        ecc = norm(e_i)

        if abs(ecc) <= 1 - 1e-6
            a = -μ / (2ξ)
        else
            throw(ArgumentError("""
                Could not convert the provided Cartesian values to Kepler elements.
                The computed eccentricity was not between 0 and 1."""
            ))
        end

        cos_i = h_i[3] / h
        cos_i = abs(cos_i) > 1 ? sign(cos_i) : cos_i
        i     = acos(cos_i)

        if abs(n) <= 1e-6
            # Equatorial orbit

            Ω = T(0)

            if abs(ecc) > 1e-6
                cos_ω = e_i[1] / ecc
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                if e_i[2] < 0
                    ω = T(2π) - ω
                end

                cos_f = dot(e_i, sr_i) / (ecc * r)
                cos_f = abs(cos_f) > 1 ? sign(cos_f) : cos_f
                f     = acos(cos_f)

                if rv < 0
                    f = T(2π) - f
                end
            else
                # Equatorial circular 

                ω = T(0)

                cos_f = sr_i[1] / r
                cos_f = abs(cos_f) > 1 ? sign(cos_f) : cos_f
                f     = acos(cos_f)

                if sr_i[2] < 0
                    f = T(2π) - f
                end
            end

        else
            # Inclined orbit
            cos_Ω = n_i[1] / n
            cos_Ω = abs(cos_Ω) > 1 ? sign(cos_Ω) : cos_Ω
            Ω     = acos(cos_Ω)

            if n_i[2] < 0
                Ω = T(2π) - Ω
            end

            if abs(ecc) < 1e-6
                # Circular

                ω = T(0)

                cos_f = dot(n_i, sr_i) / (n * r)
                cos_f = abs(cos_f) > 1 ? sign(cos_f) : cos_f
                f     = acos(cos_f)

                if sr_i[3] < 0
                    f = T(2π) - f
                end
            else

                cos_ω = dot(n_i, e_i) / (n * ecc)
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                if e_i[3] < 0
                    ω = T(2π) - ω
                end

                cos_f = dot(e_i, sr_i) / (ecc * r)
                cos_f = abs(cos_f) > 1 ? sign(cos_f) : cos_f
                f     = acos(cos_f)

                if rv < 0
                    f = T(2π) - f
                end
            end
        end
    end

    return [a, ecc, i, ω, Ω, f]
end

"""
   isapprox_angle(a, b; atol=1e-8)

   Compare two angles knowing they live in the circle
"""
isapprox_angle(a, b; atol=1e-8) = abs(rem2pi(a - b, RoundNearest)) < atol


@testset "J2 compare against poliastro" begin
    μ = 398600441800000.0
    j2 = 0.00108263
    r = 6378136.6

    cbf = TwoBodyForce(μ)
    j2f = J2Force(μ, r, j2)
    fm = ForceModel((cbf, j2f), Val(true))

    # Data from poliastro, see validation/ folder
    orbit1_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    orbit1_u1 = SA[6679.99961672, -711.28165725, -849.9863342, 1.25649055, 4.69949408, 5.93788535] .* 1e3

    orbit2_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    orbit2_u1 = SA[-5887.78600842, -1753.87582895, -2850.77957146, 3.75779874, -4.36422534, -5.06992635] .* 1e3

    orbit3_u0 = SA[63389.0685, 0, 0, 0, 1.9076657, 2.4068751] .* 1e3
    orbit3_u1 = SA[-80225.31330281, 67593.3242276, 85279.58309128, -1.64801759, -0.11879248, -0.14990463] .* 1e3

    @test propagate_orbit(fm, orbit1_u0, 90.0 * 60.0) ≈ orbit1_u1 rtol = 1e-8
    @test propagate_orbit(fm, orbit2_u0, 1440.0 * 60.0) ≈ orbit2_u1 rtol = 1e-6
    @test propagate_orbit(fm, orbit3_u0, 1440.0 * 60.0) ≈ orbit3_u1 rtol = 1e-8

end

@testset "Keplerian / Cartesian to MEE elements" begin
    @testset "Equatorial circular orbit" begin
        kepler = SA[1000.0*1e3, 0, 0, 0, 0, 0]
        mee = kepler_to_mee(kepler...)
        kepler2 = mee_to_kepler(mee...)
        @testset for i = 1:2
            @test kepler2[i] ≈ kepler[i]
        end
        @testset for i = 3:6
            @test isapprox_angle.(kepler2[i], kepler[i])
        end
    end

    @testset "Polar circular orbit" begin
        kepler = SA[1000.0*1e3, 0, π/4, 0, 0, 0]
        mee = kepler_to_mee(kepler...)
        kepler2 = mee_to_kepler(mee...)
        @testset for i = 1:2
            @test kepler2[i] ≈ kepler[i]
        end
        @testset for i = 3:6
            @test isapprox_angle.(kepler2[i], kepler[i])
        end
    end

    @testset "Random orbit" begin
        euclidean = SA[6771.358863, 1313.0, 1314.43, 0.3, 4.76807358, 6.01581168] .* 1e3
        kepler = euclid_to_kepler(euclidean..., GM_EARTH)
        mee = kepler_to_mee(kepler...)
        kepler2 = mee_to_kepler(mee...)
        mee_direct = euclid_to_mee(euclidean..., GM_EARTH)
        
        @testset "Kepler -> MEE -> Kepler" begin
            @testset for i = 1:2
                @test kepler2[i] ≈ kepler[i]
            end
            @testset for i = 3:6
                @test isapprox_angle.(kepler2[i], kepler[i])
            end
        end

        @testset "Euclidean -> MEE compared to Euclidean -> Kepler -> MEE" begin
            @testset for i = 1:5
                @test mee_direct[i] ≈ mee[i]
            end

            @test isapprox_angle.(mee_direct[6], mee[6])
        end
    end
end

@testset "Euclidean -> MEE -> Euclidean" begin
    euclidean1 = SA[6771.358863, 1313.0, 1314.43, 0.3, 4.76807358, 6.01581168] .* 1e3
    mee = euclid_to_mee(euclidean1..., GM_EARTH)
    euclidean2 = mee_to_euclid(mee..., GM_EARTH)

    @testset for i=1:6
        @test euclidean2[i] ≈ euclidean1[i]
    end
end


# The following "Edge case orbits" tests were proposed by Claude Sonnet 4.6,
# the only real "edge" case is the retrograde orbit, which is truly non-represented
# by MEE elements. They are otherwise well behaved under MEE elements!
@testset "Edge case orbits" begin
    GM = GM_EARTH

    @testset "Circular orbit (e≈0)" begin
        # ISS-like circular orbit
        r = 6771e3
        v_circ = sqrt(GM / r)
        euclidean1 = SA[r, 0.0, 0.0, 0.0, v_circ, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] atol = 1e-6
        end
    end

    @testset "Polar orbit (i=90°)" begin
        euclidean1 = SA[7000e3, 0.0, 0.0, 0.0, 0.0, sqrt(GM / 7000e3)]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] atol = 1e-6
        end
    end

    @testset "Equatorial prograde orbit (i≈0)" begin
        r = 42164e3  # GEO radius
        v = sqrt(GM / r)
        euclidean1 = SA[r, 0.0, 0.0, 0.0, v, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] atol = 1e-6
        end
    end

    @testset "High eccentricity elliptical (e≈0.9, Molniya-like)" begin
        # Molniya orbit: high apogee, low perigee
        rp = 6778e3        # perigee ~400 km
        ra = 40000e3       # apogee
        a  = (rp + ra) / 2
        e  = (ra - rp) / (ra + rp)
        vp = sqrt(GM * (2/rp - 1/a))
        euclidean1 = SA[rp, 0.0, 0.0, 0.0, vp * cosd(63.4), vp * sind(63.4)]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] rtol = 1e-8
        end
    end

    @testset "Parabolic trajectory (e=1)" begin
        r = 7000e3
        v_escape = sqrt(2GM / r)  # exactly escape velocity → e = 1
        euclidean1 = SA[r, 0.0, 0.0, 0.0, v_escape, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] rtol = 1e-8
        end
    end

    @testset "Hyperbolic trajectory (e>1, flyby)" begin
        r = 7000e3
        v_hyp = 1.3 * sqrt(2GM / r)  # excess velocity beyond escape
        euclidean1 = SA[r, 0.0, 0.0, 0.0, v_hyp, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] rtol = 1e-8
        end
    end

    @testset "Near-degenerate rectilinear (very low angular momentum)" begin
        # Almost radial trajectory — tiny tangential velocity
        r = 7000e3
        v_circ = sqrt(GM / r)
        euclidean1 = SA[r, 0.0, 0.0, 1.0, v_circ * 0.001, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] rtol = 1e-6
        end
    end

    @testset "Sun-synchronous-like (i≈97.8°)" begin
        euclidean1 = SA[6978e3, 1313.0, 1314.43, 0.3, 4.76807358, 6.01581168] .* 1e0
        # Patch inclination by rotating velocity into ~97.8° plane
        inc = deg2rad(97.8)
        r   = 6978e3
        v   = sqrt(GM / r)
        euclidean1 = SA[r, 0.0, 0.0, 0.0, v * cos(inc), v * sin(inc)]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test euclidean2[i] ≈ euclidean1[i] rtol = 1e-8
        end
    end

    @testset "Retrograde equatorial orbit (i=180°) is not representable in MEE" begin
        r = 7500e3
        v = sqrt(GM / r)
        euclidean1 = SA[r, 0.0, 0.0, 0.0, -v, 0.0]
        mee = euclid_to_mee(euclidean1..., GM)
        euclidean2 = mee_to_euclid(mee..., GM)
        @testset for i = 1:6
            @test isnan(euclidean2[i])
        end
    end

end

@testset "Twobody kepler vs newton" begin
    orbit_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3

    # NEWTON propagation
    orbit_u1 = propagate_orbit(EARTH_FM_NEWTON, orbit_u0, 3600.0)

    # KEPLER propagation
    orbit_u0_kepler = euclid_to_kepler(orbit_u0..., GM_EARTH)
    orbit_u0_mee = kepler_to_mee(orbit_u0_kepler...)

    orbit_u1_mee = propagate_orbit(EARTH_FM_KEPLER, orbit_u0_mee, 3600.0)
    orbit_u1_mee_kepler = mee_to_kepler(orbit_u1_mee...)
    orbit_u1_mee_euclid = kepler_to_euclid(orbit_u1_mee_kepler..., GM_EARTH)

    @test orbit_u1 ≈ orbit_u1_mee_euclid rtol=1e-6
end


@testset "J2 kepler vs newton" begin
    orbit_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3

    # NEWTON propagation
    orbit_u1 = propagate_orbit(EARTH_FM_WITH_J2_NEWTON, orbit_u0, 3600.0)

    # KEPLER propagation
    orbit_u0_kepler = euclid_to_kepler(orbit_u0..., GM_EARTH)
    orbit_u0_mee = kepler_to_mee(orbit_u0_kepler...)

    orbit_u1_mee = propagate_orbit(EARTH_FM_WITH_J2_KEPLER, orbit_u0_mee, 3600.0)
    orbit_u1_mee_kepler = mee_to_kepler(orbit_u1_mee...)
    orbit_u1_mee_euclid = kepler_to_euclid(orbit_u1_mee_kepler..., GM_EARTH)

    @test orbit_u1 ≈ orbit_u1_mee_euclid rtol=1e-6
end


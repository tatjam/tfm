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

@testset "Keplerian to MEE elements" begin
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
        @info kepler
        @info kepler2
        @testset for i = 1:2
            @test kepler2[i] ≈ kepler[i]
        end
        @testset for i = 3:6
            @test isapprox_angle.(kepler2[i], kepler[i])
        end
    end

    @testset "Random orbit" begin
        euclidean = SA[6771.358863, 1313.0, 1314.43, 0.3, 4.76807358, 6.01581168] .* 1e3
        kepler = kepler_to_array(rv_to_kepler(euclidean[1:3], euclidean[4:6], GM_EARTH))
        mee = kepler_to_mee(kepler...)
        kepler2 = mee_to_kepler(mee...)
        @info kepler
        @info kepler2
        @testset for i = 1:2
            @test kepler2[i] ≈ kepler[i]
        end
        @testset for i = 3:6
            @test isapprox_angle.(kepler2[i], kepler[i])
        end
    end
end

@testset "J2 kepler vs newton" begin
    orbit1_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3

end


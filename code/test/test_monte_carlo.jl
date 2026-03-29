@testset "Monte carlo of various samples = individual propagation of various samples" begin
    samples = [
        SA[1000.0, 0.0, 0.0, 1000.0, 0.0, 0.0],
        SA[5435.0, 5435.0, 4242.0, 4324.0, 4324.0, 4244.0],
        SA[1555.0, 5555.0, 5555.0, 1234.0, 1234.0, 0.0],
    ]

    results = run_monte_carlo(EARTH_FM_WITH_J2_NEWTON, samples, 1000)

    for i in eachindex(samples)
        ind_result = propagate_orbit(EARTH_FM_WITH_J2_NEWTON, samples[i], 1000)
        @test ind_result ≈ results[i] rtol = 1e-14
    end
end
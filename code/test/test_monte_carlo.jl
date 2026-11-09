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

@testset "Monte carlo with extras (log β)" begin
    drag = SimpleDragForce(1e-12, 6771e3, 8e4)
    fm = ForceModel(GM_EARTH, (drag,), Val(true))
    logβ = log(50.0)
    samples = [
        SA[4789e3, 0.0, 4789e3, 0.0, 7.5e3, 0.0, logβ],
        SA[6771e3, 0.0, 0.0, 0.0, 4.768e3, 6.016e3, logβ],
    ]

    results = run_monte_carlo(fm, samples, 1000)

    for i in eachindex(samples)
        ind_result = propagate_orbit(fm, samples[i], 1000)
        @test ind_result ≈ results[i] rtol = 1e-14
    end
end

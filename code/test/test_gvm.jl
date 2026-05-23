@testset "Propagation" begin

    @testset "Very short orbital case matches sufficient Monte Carlo samples" begin
        Δt = 60
        nsamples = 10000

        orbit_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
        orbit_u0_kepler = euclid_to_kepler(orbit_u0..., GM_EARTH)
        orbit_u0_mee = kepler_to_mee(orbit_u0_kepler...)

        P = I(5) * 0.001
        β = [0.0, 0.0, 0.0, 0.0, 0.0]
        Γ = zeros(5, 5)
        κ = 0.1
        
        dist = GaussVonMises(orbit_u0_mee[1:5], orbit_u0_mee[6], β, Γ, κ, P=P)
        # A sufficiently short Keplerian problem should be very exactly modelled
        end_dist = run_gvm(EARTH_FM_WITH_J2_KEPLER, dist, Δt)

        # Monte-Carlo samples
        samples = rand(dist, nsamples)
        samples_mc = run_monte_carlo(EARTH_FM_WITH_J2_KEPLER, samples, Δt)

        # Flatten to a conventional matrix 
        samples_mc_mat = reduce(hcat, samples_mc)

        # TODO: Fit the resulting MC samples to a GVM
        # TODO: Expect them to reasonably match

    end

end

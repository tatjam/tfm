# EnergyTestOverPropagation.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox

function run_comparison(fm, fm_kepler, starting_dist, starting_dist_kep, t1, Δt)
    # Samples of ground-truth Monte Carlo
    NSAMPLES_TRUTH = 10000
    # Samples of methods evaluated against ground-truth
    NSAMPLES = 1000

    # We propagate a 10000 point Monte Carlo as "ground-truth" distribution
    ref_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]

    # This is a lighter Mont Carlo just to compare it to the other methods
    mc_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, NSAMPLES))]

    current_ref_samples = copy(ref_samples)
    current_mc_samples = copy(mc_samples)

    μ₀ = mean(starting_dist)
    P₀ = cov(starting_dist)

    μ₀_kepler = mean(starting_dist_kep)
    P₀_kepler = cov(starting_dist_kep)

    nsteps = floor(t1 / Δt)

    # mc, ut, stm, ut_kep, stm_kep
    energies = zeros(Float32, (5, nsteps))

    for t = 0:Δt:t1
        t = Δt * i

        # Propagate Monte Carlo from t - Δt -> t
        # WARNING: There's a clear issue here, if the force model is not time-independent
        # we need to update its internal time for the Monte Carlo step, as it doesn't
        # start at t=0 but at t=t implicitly!
        current_ref_samples = run_monte_carlo(fm, current_ref_samples, Δt)
        current_mc_samples = run_monte_carlo(fm, current_mc_samples, Δt)

        # Propagate statistical propagators from 0 -> t
        ut_dist = run_ut(fm, μ₀, P₀, t)
        stm_dist = run_stm(fm, μ₀, P₀, t)
        ut_dist_kepler = run_ut(fm_kepler, μ₀_kepler, P₀_kepler, t)
        stm_dist_kepler = run_stm(fm_kepler, μ₀_kepler, P₀_kepler, t)

        # Sample (and transform to Euclidean coords if needed)
        ut_samples = stack(rand(ut_dist, NSAMPLES))
        stm_samples = stack(rand(stm_dist, NSAMPLES))
        ut_samples_kepler = mapslices(
            v -> mee_to_euclid.(v..., GM_EARTH),
            stack(rand(ut_dist_kepler, NSAMPLES)),
            dims = 1,
        )
        stm_samples_kepler = mapslices(
            v -> mee_to_euclid.(v..., GM_EARTH),
            stack(rand(stm_dist_kepler, NSAMPLES)),
            dims = 1,
        )

        @info "Propagation to t = ", t, " complete."

        # TODO: Whiten?
        energy_mc = energy_metric(current_mc_samples, current_ref_samples, Float32)
        energy_ut = energy_metric(ut_samples, current_ref_samples, Float32)
        energy_stm = energy_metric(stm_samples, current_ref_samples, Float32)
        energy_ut_kepler = energy_metric(ut_samples_kepler, current_ref_samples, Float32)
        energy_stm_kepler = energy_metric(stm_samples_kepler, current_ref_samples, Float32)

        @info "Energies of t = ", t, " complete."

        energies[0, i] = energy_mc
        energies[1, i] = energy_ut
        energies[2, i] = energy_stm
        energies[3, i] = energy_ut_kepler
        energies[4, i] = energy_stm_kepler
    end

    # TODO: Save data for posterior plotting
end



function main()

    END_T = 3600.0 * 12.0
    DELTA_T = 60.0

    fm = EARTH_FM_WITH_J2_NEWTON
    fmk = EARTH_FM_WITH_J2_KEPLER
    μ = [9000e3, 0, 0, 0, 6620, 0]
    σ = Diagonal([100e3, 10e3, 10e3, 1.0, 100.0, 1.0])

    starting_dist = MvNormal(μ, σ^2)

    # Transform to MEE coordinates, assuming normal after the non-linear transform
    starting_dist_kep = ut_propagate(v -> euclid_to_mee(v..., GM_EARTH), μ, σ^2, α = 1e-1)
    run_comparison(fm, fmk, starting_dist, starting_dist_kep, END_T, DELTA_T)

end

main()


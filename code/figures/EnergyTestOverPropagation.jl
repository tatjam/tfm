# EnergyTestOverPropagation.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox
using HDF5


function whitener(ref)
    ref_mean = mean(ref, dims = 2)
    ref_cov = cov(ref, dims = 2)
    L = cholesky(Symmetric(ref_cov)).L

    return X -> L \ (X .- ref_mean)
end

function run_comparison(
    fm,
    fm_kepler,
    starting_dist,
    starting_dist_kep,
    t1,
    Δt,
    outfile = "comparison.h5",
)
    # Samples of ground-truth Monte Carlo
    NSAMPLES_TRUTH = 3000
    # Samples of methods evaluated against ground-truth
    NSAMPLES = 800
    NMETHODS = 6

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

    ts = 0:Δt:t1
    nsteps = length(ts)

    h5open(outfile, "w") do f

        attrs(f)["NSAMPLES"] = NSAMPLES
        attrs(f)["NSAMPLES_TRUTH"] = NSAMPLES_TRUTH
        attrs(f)["delta_t"] = Δt
        attrs(f)["methods"] = "mc,ut,stm,ut_kep,stm_kep,ut_kep_alt"
        f["ts"] = collect(ts)

        d_truth = create_dataset(
            f,
            "truth",
            Float32,
            (6, NSAMPLES_TRUTH, nsteps);
            chunk = (6, NSAMPLES_TRUTH, 1),
            deflate = 3,
        )

        d_samples = create_dataset(
            f,
            "samples",
            Float32,
            (6, NSAMPLES, NMETHODS, nsteps);
            chunk = (6, NSAMPLES, NMETHODS, 1),
            deflate = 3,
        )

        d_energy = create_dataset(f, "energies", Float32, (NMETHODS, nsteps))

        for i = 1:nsteps
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
            # The MEE covariance mixes units spanning ~16 orders of magnitude,
            # it must be equilibrated before any eigendecomposition or Cholesky
            ut_dist_kepler = run_ut(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true)
            stm_dist_kepler =
                run_stm(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true)

            # UT but with sigma points generated in the euclidean space (experiment)
            sigma_eucl = SigmaVectors(μ₀, P₀, 1e-3, 0, 2)
            sigma_mee = SigmaVectors(
                SMatrix(
                    mapslices(v -> euclid_to_mee(v..., GM_EARTH), sigma_eucl.χ, dims = 1),
                ),
                sigma_eucl.W,
                sigma_eucl.W0c,
            )
            ut_dist_kepler_alt =
                ut_propagate_sigma(v -> propagate_orbit(fm_kepler, v, t), sigma_mee)

            # Sample (and transform to Euclidean coords if needed)
            ut_samples = stack(rand(ut_dist, NSAMPLES))
            stm_samples = stack(rand(stm_dist, NSAMPLES))
            ut_samples_kepler = mapslices(
                v -> mee_to_euclid(v..., GM_EARTH),
                stack(rand(ut_dist_kepler, NSAMPLES)),
                dims = 1,
            )
            stm_samples_kepler = mapslices(
                v -> mee_to_euclid(v..., GM_EARTH),
                stack(rand(stm_dist_kepler, NSAMPLES)),
                dims = 1,
            )
            ut_samples_kepler_alt = mapslices(
                v -> mee_to_euclid(v..., GM_EARTH),
                stack(rand(ut_dist_kepler_alt, NSAMPLES)),
                dims = 1,
            )

            @info "Propagation complete" t

            mc_mat = reduce(hcat, current_mc_samples)
            ref_mat = reduce(hcat, current_ref_samples)

            d_truth[:, :, i] = ref_mat
            d_samples[:, :, 1, i] = mc_mat
            d_samples[:, :, 2, i] = ut_samples
            d_samples[:, :, 3, i] = stm_samples
            d_samples[:, :, 4, i] = ut_samples_kepler
            d_samples[:, :, 5, i] = stm_samples_kepler
            d_samples[:, :, 6, i] = ut_samples_kepler_alt

            # Withen with the reference
            w = whitener(ref_mat)
            ref_white = w(ref_mat)

            energy = X -> energy_metric(w(X), ref_white; precision = Float32)

            energy_mc = energy(mc_mat)
            energy_ut = energy(ut_samples)
            energy_stm = energy(stm_samples)
            energy_ut_kepler = energy(ut_samples_kepler)
            energy_stm_kepler = energy(stm_samples_kepler)
            energy_ut_kepler_alt = energy(ut_samples_kepler_alt)

            @info "Energies complete" t

            d_energy[:, i] = Float32[
                energy_mc,
                energy_ut,
                energy_stm,
                energy_ut_kepler,
                energy_stm_kepler,
                energy_ut_kepler_alt,
            ]

            flush(f)
        end
    end

    return outfile
end

function main()
    END_T = 3600.0 * 24.0 * 10
    DELTA_T = 6000.0

    fm = EARTH_FM_WITH_J2_NEWTON
    fmk = EARTH_FM_WITH_J2_KEPLER
    μ = [9000e3, 0, 0, 0, 8620, 0]
    σ = Diagonal([1e3, 1e3, 1e3, 0.02, 0.02, 0.02])

    starting_dist = MvNormal(μ, σ^2)

    # Transform to MEE coordinates, assuming normal after the non-linear transform
    # (equilibrate: the output MEE covariance mixes units, see run_ut/ut_propagate)
    starting_dist_kep = ut_propagate(
        v -> euclid_to_mee(v..., GM_EARTH),
        μ,
        σ^2,
        α = 1e-1,
        equilibrate = true,
    )

    ts, truth, sampled_states, energies =
        run_comparison(fm, fmk, starting_dist, starting_dist_kep, END_T, DELTA_T)
end

main()

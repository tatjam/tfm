# MahalanobisTest.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using Distributions
using Statistics
using StaticArrays
using LinearAlgebra
using SatelliteToolbox
using HDF5

const NSAMPLES_TRUTH = 10000
const NSAMPLES_MC_MAPPING = 1000
const METHODS = [
    (name = "UT", color = :dodgerblue3),
    (name = "STM", color = :chocolate1),
    (name = "UT (MEE)", color = :deepskyblue),
    (name = "STM (MEE)", color = :crimson),
    (name = "UT (MEE->EUCL UT)", color = :teal),
    (name = "STM (MEE->EUCL UT)", color = :darkorchid),
    (name = "UT (MEE->EUCL MC)", color = :cadetblue),
    (name = "STM (MEE->EUCL MC)", color = :mediumvioletred),
]
const QUANTILES = [0.05, 0.5, 0.95]

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

    # Ground truth distribution
    ref_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]

    current_ref_samples = copy(ref_samples)

    μ₀ = mean(starting_dist)
    P₀ = cov(starting_dist)

    μ₀_kepler = mean(starting_dist_kep)
    P₀_kepler = cov(starting_dist_kep)

    ts = 0:Δt:t1
    nsteps = length(ts)

    h5open(outfile, "w") do f

        attrs(f)["NSAMPLES_TRUTH"] = NSAMPLES_TRUTH
        attrs(f)["delta_t"] = Δt
        attrs(f)["methods"] = map(v->v.name, METHODS)
        attrs(f)["QUANTILES"] = QUANTILES
        f["ts"] = collect(ts)

        d_truth = create_dataset(
            f,
            "truth",
            Float32,
            (6, NSAMPLES_TRUTH, nsteps);
            chunk = (6, NSAMPLES_TRUTH, 1),
            deflate = 3,
        )


        d_means = create_dataset(
            f,
            "means",
            Float32,
            (6, length(METHODS), nsteps);
            chunk = (6, length(METHODS), 1),
            deflate = 3,
        )

        d_covsm = create_dataset(
            f,
            "covs",
            Float32,
            (6, 6, length(METHODS), nsteps);
            chunk = (6, 6, length(METHODS), 1),
            deflate = 3,
        )

        d_mahalanobis_tot =
            create_dataset(f, "mahalanobis_tot", Float32, (length(METHODS), nsteps))
        d_mahalanobis_q = create_dataset(
            f,
            "mahalanobis_q",
            Float32,
            (length(QUANTILES), length(METHODS), nsteps),
        )

        for i = 1:nsteps
            t = Δt * i
            current_ref_samples = run_monte_carlo(fm, current_ref_samples, Δt)
            current_ref_samples_kepler =
                map(v -> euclid_to_mee(v..., GM_EARTH), current_ref_samples)

            mc_mat = reduce(hcat, current_ref_samples)
            mc_mat_kepler = reduce(hcat, current_ref_samples_kepler)

            d_truth[:, :, i] = mc_mat

            ut_dist = run_ut(fm, μ₀, P₀, t)
            stm_dist = run_stm(fm, μ₀, P₀, t)
            ut_dist_kepler = run_ut(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true, α=1e-1)
            stm_dist_kepler =
                run_stm(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true)

            d_means[:, 1, i] = mean(ut_dist)
            d_covsm[:, :, 1, i] = cov(ut_dist)

            d_means[:, 2, i] = mean(stm_dist)
            d_covsm[:, :, 2, i] = cov(stm_dist)

            d_means[:, 3, i] = mean(ut_dist_kepler)
            d_covsm[:, :, 3, i] = cov(ut_dist_kepler)

            d_means[:, 4, i] = mean(stm_dist_kepler)
            d_covsm[:, :, 4, i] = cov(stm_dist_kepler)

            # MEE -> Euclid mapping via UT
            ut_dist_kepler_eucl_ut = ut_propagate(
                v -> mee_to_euclid(v..., GM_EARTH),
                mean(ut_dist_kepler),
                cov(ut_dist_kepler),
                α = 1e-1,
                equilibrate = true,
            )
            stm_dist_kepler_eucl_ut = ut_propagate(
                v -> mee_to_euclid(v..., GM_EARTH),
                mean(stm_dist_kepler),
                cov(stm_dist_kepler),
                α = 1e-1,
                equilibrate = true,
            )

            d_means[:, 5, i] = mean(ut_dist_kepler_eucl_ut)
            d_covsm[:, :, 5, i] = cov(ut_dist_kepler_eucl_ut)

            d_means[:, 6, i] = mean(stm_dist_kepler_eucl_ut)
            d_covsm[:, :, 6, i] = cov(stm_dist_kepler_eucl_ut)

            # MEE -> Euclid mapping via sampling
            samples_ut_dist_kepler = rand(ut_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_ut_dist_kepler_eucl = mapslices(
                v->mee_to_euclid(v..., GM_EARTH),
                samples_ut_dist_kepler,
                dims = 1,
            )
            ut_dist_kepler_eucl_mc = fit(MvNormal, samples_ut_dist_kepler_eucl)

            d_means[:, 7, i] = mean(ut_dist_kepler_eucl_mc)
            d_covsm[:, :, 7, i] = cov(ut_dist_kepler_eucl_mc)


            samples_stm_dist_kepler = rand(stm_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_stm_dist_kepler_eucl = mapslices(
                v->mee_to_euclid(v..., GM_EARTH),
                samples_stm_dist_kepler,
                dims = 1,
            )
            stm_dist_kepler_eucl_mc = fit(MvNormal, samples_stm_dist_kepler_eucl)
            d_means[:, 8, i] = mean(stm_dist_kepler_eucl_mc)
            d_covsm[:, :, 8, i] = cov(stm_dist_kepler_eucl_mc)

            # Mahalanobis distance of MC samples to Euclidean methods
            mahal_ut = sqmahal(ut_dist, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[1, i] = sum(mahal_ut)
            d_mahalanobis_q[:, 1, i] = quantile(mahal_ut, QUANTILES)

            mahal_stm = sqmahal(stm_dist, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[2, i] = sum(mahal_stm)
            d_mahalanobis_q[:, 2, i] = quantile(mahal_stm, QUANTILES)

            mahal_ut_kep_eucl_ut = sqmahal(ut_dist_kepler_eucl_ut, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[5, i] = sum(mahal_ut_kep_eucl_ut)
            d_mahalanobis_q[:, 5, i] = quantile(mahal_ut_kep_eucl_ut, QUANTILES)

            mahal_stm_kep_eucl_ut =
                sqmahal(stm_dist_kepler_eucl_ut, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[6, i] = sum(mahal_stm_kep_eucl_ut)
            d_mahalanobis_q[:, 6, i] = quantile(mahal_stm_kep_eucl_ut, QUANTILES)

            mahal_ut_kep_eucl_mc = sqmahal(ut_dist_kepler_eucl_mc, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[7, i] = sum(mahal_ut_kep_eucl_mc)
            d_mahalanobis_q[:, 7, i] = quantile(mahal_ut_kep_eucl_mc, QUANTILES)

            mahal_stm_kep_eucl_mc =
                sqmahal(stm_dist_kepler_eucl_mc, mc_mat) / NSAMPLES_TRUTH
            d_mahalanobis_tot[8, i] = sum(mahal_stm_kep_eucl_mc)
            d_mahalanobis_q[:, 8, i] = quantile(mahal_stm_kep_eucl_mc, QUANTILES)

            # Fit the Von Mises distribution to the wrapped normal on MEE methods
            L_COORD = 6
            ut_kep_mahal_gvm = moment_match(ut_dist_kepler, L_COORD)
            stm_kep_mahal_gvm = moment_match(stm_dist_kepler, L_COORD)

            # Von-Mises Mahalanobis distance of MC samples to MEE methods
            mahal_ut_kep_gvm =
                mapslices(v->mahalanobis(v, ut_kep_mahal_gvm), mc_mat_kepler, dims = 1) /
                NSAMPLES_TRUTH
            d_mahalanobis_tot[3, i] = sum(mahal_ut_kep_gvm)
            d_mahalanobis_q[:, 3, i] = quantile(mahal_ut_kep_gvm, QUANTILES)

            mahal_stm_kep_gvm =
                mapslices(v->mahalanobis(v, stm_kep_mahal_gvm), mc_mat_kepler, dims = 1) /
                NSAMPLES_TRUTH
            d_mahalanobis_tot[4, i] = sum(mahal_stm_kep_gvm)
            d_mahalanobis_q[:, 4, i] = quantile(mahal_stm_kep_gvm, QUANTILES)

            @info "Propagation complete" t
            flush(f)
        end


    end

    return outfile
end

function main()
    END_T = 3600.0 * 24.0 * 2
    DELTA_T = 360.0

    fm = EARTH_FM_WITH_J2_NEWTON
    fmk = EARTH_FM_WITH_J2_KEPLER
    μ = [9000e3, 0, 0, 0, 8620, 0]
    σ = Diagonal([1e3, 1e3, 1e3, 0.02, 0.02, 0.02])

    starting_dist = MvNormal(μ, σ^2)

    # Transform to MEE coordinates, assuming normal after the non-linear transform
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

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
    (name = "Ref MC", color = :black),
    (name = "Resampled MC", color = :gray),
]
const QUANTILES = [0.05, 0.5, 0.95]

const METHOD_IDX = Dict(m.name => j for (j, m) in enumerate(METHODS))

function save_dist!(d_mean, d_covm, name, i, dist)
    j = METHOD_IDX[name]
    d_mean[:, j, i] = mean(dist)
    d_covm[:, :, j, i] = cov(dist)
end

function save_mahal!(d_q, name, i, dist, mc_mat)
    d = sqmahal(dist, mc_mat)
    j = METHOD_IDX[name]
    d_q[:, j, i] = quantile(vec(d), QUANTILES)
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

    # Ground truth distribution
    ref_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]
    current_ref_samples = copy(ref_samples)

    # Resampling to have a baseline
    alt_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]
    current_alt_samples = copy(alt_samples)

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

        d_mahal_q = create_dataset(
            f,
            "mahalanobis_q",
            Float32,
            (length(QUANTILES), length(METHODS), nsteps),
        )

        for i = 1:nsteps
            t = Δt * i
            current_ref_samples = run_monte_carlo(fm, current_ref_samples, Δt)
            current_alt_samples = run_monte_carlo(fm, current_alt_samples, Δt)
            current_ref_samples_kepler =
                map(v -> euclid_to_mee(v..., GM_EARTH), current_ref_samples)

            mc_mat = reduce(hcat, current_ref_samples)
            mc_alt_mat = reduce(hcat, current_alt_samples)
            mc_mat_kepler = reduce(hcat, current_ref_samples_kepler)

            d_truth[:, :, i] = mc_mat

            ut_dist = run_ut(fm, μ₀, P₀, t; equilibrate=true)
            stm_dist = run_stm(fm, μ₀, P₀, t; equilibrate=true)
            ut_dist_kepler =
                run_ut(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true, α = 1e-1)
            stm_dist_kepler =
                run_stm(fm_kepler, μ₀_kepler, P₀_kepler, t; equilibrate = true)

            save_dist!(d_means, d_covsm, "UT", i, ut_dist)
            save_dist!(d_means, d_covsm, "STM", i, stm_dist)

            save_dist!(d_means, d_covsm, "UT (MEE)", i, ut_dist_kepler)
            save_dist!(d_means, d_covsm, "STM (MEE)", i, stm_dist_kepler)
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

            save_dist!(d_means, d_covsm, "UT (MEE->EUCL UT)", i, ut_dist_kepler_eucl_ut)
            save_dist!(d_means, d_covsm, "STM (MEE->EUCL UT)", i, stm_dist_kepler_eucl_ut)

            # MEE -> Euclid mapping via sampling
            samples_ut_dist_kepler = rand(ut_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_ut_dist_kepler_eucl = mapslices(
                v->mee_to_euclid(v..., GM_EARTH),
                samples_ut_dist_kepler,
                dims = 1,
            )
            ut_dist_kepler_eucl_mc = fit(MvNormal, samples_ut_dist_kepler_eucl)

            save_dist!(d_means, d_covsm, "UT (MEE->EUCL MC)", i, ut_dist_kepler_eucl_mc)

            samples_stm_dist_kepler = rand(stm_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_stm_dist_kepler_eucl = mapslices(
                v->mee_to_euclid(v..., GM_EARTH),
                samples_stm_dist_kepler,
                dims = 1,
            )
            stm_dist_kepler_eucl_mc = fit(MvNormal, samples_stm_dist_kepler_eucl)
            save_dist!(d_means, d_covsm, "STM (MEE->EUCL MC)", i, stm_dist_kepler_eucl_mc)

            # Reference and alternate sampling MC
            ref_mc = fit(MvNormal, mc_mat)
            alt_mc = fit(MvNormal, mc_alt_mat)

            save_dist!(d_means, d_covsm, "Ref MC", i, ref_mc)
            save_dist!(d_means, d_covsm, "Resampled MC", i, alt_mc)

            # Mahalanobis distance of MC samples to Euclidean methods
            save_mahal!(d_mahal_q, "UT", i, ut_dist, mc_mat)
            save_mahal!(d_mahal_q, "STM", i, stm_dist, mc_mat)
            save_mahal!(
                d_mahal_q,
                "UT (MEE->EUCL UT)",
                i,
                ut_dist_kepler_eucl_ut,
                mc_mat,
            )
            save_mahal!(
                d_mahal_q,
                "STM (MEE->EUCL UT)",
                i,
                stm_dist_kepler_eucl_ut,
                mc_mat,
            )
            save_mahal!(
                d_mahal_q,
                "UT (MEE->EUCL MC)",
                i,
                ut_dist_kepler_eucl_mc,
                mc_mat,
            )
            save_mahal!(
                d_mahal_q,
                "STM (MEE->EUCL MC)",
                i,
                stm_dist_kepler_eucl_mc,
                mc_mat,
            )

            save_mahal!(d_mahal_q, "Ref MC", i, ref_mc, mc_mat)
            save_mahal!(d_mahal_q, "Resampled MC", i, alt_mc, mc_mat)

            # Fit the Von Mises distribution to the wrapped normal on MEE methods
            L_COORD = 6
            ut_kep_mahal_gvm = moment_match(ut_dist_kepler, L_COORD)
            stm_kep_mahal_gvm = moment_match(stm_dist_kepler, L_COORD)

            # Von-Mises Mahalanobis distance of MC samples to MEE methods
            mahal_ut_kep_gvm =
                mapslices(v->mahalanobis(v, ut_kep_mahal_gvm), mc_mat_kepler, dims = 1)
            d_mahal_q[:, 3, i] = quantile(mahal_ut_kep_gvm, QUANTILES)

            mahal_stm_kep_gvm =
                mapslices(v->mahalanobis(v, stm_kep_mahal_gvm), mc_mat_kepler, dims = 1)
            d_mahal_q[:, 4, i] = quantile(mahal_stm_kep_gvm, QUANTILES)

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

    run_comparison(fm, fmk, starting_dist, starting_dist_kep, END_T, DELTA_T)
end

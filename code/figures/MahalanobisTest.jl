# MahalanobisTest.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using Distributions
using Statistics
using StaticArrays
using LinearAlgebra
using Random
using SatelliteToolbox
using HDF5

const NSAMPLES_TRUTH = 10000
const NSAMPLES_MC_MAPPING = 100_000

# nothing = disable drag
const DRAG_LOGβ_MEAN = nothing
# const DRAG_LOGβ_MEAN = log(5.0)
const DRAG_LOGβ_VAR = 0.09
const DRAG_ρ₀ = 1e-10  # kg/m³
const DRAG_r₀ = 6771e3 # m
const DRAG_H = 8e4     # m

const L_STATE = DRAG_LOGβ_MEAN === nothing ? 6 : 7
# The angular coordinate is the last element of the state either way
const L_COORD = L_STATE

const METHODS = [
    (name = "UT", color = :dodgerblue3),
    (name = "STM", color = :chocolate1),
    (name = "UT (MEE)", color = :deepskyblue),
    (name = "STM (MEE)", color = :crimson),
    (name = "UT (MEE->EUCL UT)", color = :teal),
    (name = "STM (MEE->EUCL UT)", color = :darkorchid),
    (name = "UT (MEE->EUCL MC)", color = :cadetblue),
    (name = "STM (MEE->EUCL MC)", color = :mediumvioletred),
    (name = "GVM :ad (MEE->EUCL MC)", color = :forestgreen),
    (name = "GVM :nodes (MEE->EUCL MC)", color = :olivedrab),
    (name = "GVM :ad refine (MEE->EUCL MC)", color = :darkgreen),
    (name = "GVM :nodes refine (MEE->EUCL MC)", color = :limegreen),
    (name = "Ref MC", color = :black),
    (name = "Resampled MC", color = :gray),
    (name = "UT moment-match GVM (MEE->EUCL MC)", color = :goldenrod),
]
const QUANTILES = [0.05, 0.5, 0.95]

const METHOD_IDX = Dict(m.name => j for (j, m) in enumerate(METHODS))

const GVM_VARIANTS = [
    (name = "GVM :ad (MEE->EUCL MC)", mee_name = "GVM :ad", slopes = :ad, refine = false),
    (
        name = "GVM :nodes (MEE->EUCL MC)",
        mee_name = "GVM :nodes",
        slopes = :nodes,
        refine = false,
    ),
    (
        name = "GVM :ad refine (MEE->EUCL MC)",
        mee_name = "GVM :ad refine",
        slopes = :ad,
        refine = true,
    ),
    (
        name = "GVM :nodes refine (MEE->EUCL MC)",
        mee_name = "GVM :nodes refine",
        slopes = :nodes,
        refine = true,
    ),
]

const MEE_METHODS = [
    "UT (MEE)",
    "STM (MEE)",
    "UT (MEE) GVM moment-match",
    "STM (MEE) GVM moment-match",
    "GVM :ad",
    "GVM :nodes",
    "GVM :ad refine",
    "GVM :nodes refine",
]
const MEE_IDX = Dict(n => j for (j, n) in enumerate(MEE_METHODS))

# Appends the (independent) log β dimension to a distribution when drag is on.
# For MEE-space distributions (mee = true) the extra is inserted *before* the
# angular element, matching the [p, f, g, h, k, extras..., L] state layout.
function with_drag(dist::MvNormal; mee = false)
    if DRAG_LOGβ_MEAN === nothing
        return dist
    end
    μ = mean(dist)
    P = Matrix(cov(dist))
    n = length(μ)
    if !mee
        return MvNormal(vcat(μ, DRAG_LOGβ_MEAN),
                        cat(P, fill(DRAG_LOGβ_VAR, 1, 1); dims = (1, 2)))
    end
    # Old indices 1:n-1 keep their place, old n (angular) shifts to n+1,
    # and the new index n is the log β extra
    μ_new = vcat(μ[1:(end-1)], DRAG_LOGβ_MEAN, μ[end])
    P_new = zeros(n + 1, n + 1)
    P_new[1:(n-1), 1:(n-1)] .= P[1:(n-1), 1:(n-1)]
    P_new[1:(n-1), n+1]     .= P[1:(n-1), n]
    P_new[n+1, 1:(n-1)]     .= P[n, 1:(n-1)]
    P_new[n+1, n+1]         = P[n, n]
    P_new[n, n]             = DRAG_LOGβ_VAR
    return MvNormal(μ_new, P_new)
end

# Sample conversions between Newton and MEE coordinates, carrying the log β
# extra, which sits before the angular element in MEE coordinates
function newton_to_mee(v)
    if DRAG_LOGβ_MEAN === nothing
        return euclid_to_mee(v..., GM_EARTH)
    end
    mee = euclid_to_mee(v[SA[1, 2, 3, 4, 5, 6]]..., GM_EARTH)
    return SA[mee[SA[1, 2, 3, 4, 5]]..., v[end], mee[6]]
end
mee_to_newton(v) = DRAG_LOGβ_MEAN === nothing ? mee_to_euclid(v..., GM_EARTH) :
    SA[mee_to_euclid(v[SA[1, 2, 3, 4, 5]]..., v[end], GM_EARTH)..., v[6]]

function save_dist!(d_mean, d_covm, name, i, dist)
    j = METHOD_IDX[name]
    d_mean[:, j, i] = mean(dist)
    d_covm[:, :, j, i] = cov(dist)
end

function save_nan!(d_mean, d_covm, d_q, name, i)
    j = METHOD_IDX[name]
    d_mean[:, j, i] = fill(Float32(NaN), L_STATE)
    d_covm[:, :, j, i] = fill(Float32(NaN), L_STATE, L_STATE)
    d_q[:, j, i] = fill(Float32(NaN), length(QUANTILES))
end

function save_mahal!(d_q, name, i, dist, mc_mat)
    d = sqmahal(dist, mc_mat)
    j = METHOD_IDX[name]
    d_q[:, j, i] = quantile(vec(d), QUANTILES)
end

unwrap_to(L, ref) = ref + rem2pi(L - ref, RoundNearest)

function mahal_mee(dist::MvNormal, mat)
    m = Matrix{Float64}(mat)
    m[L_COORD, :] .= unwrap_to.(m[L_COORD, :], mean(dist)[L_COORD])
    return vec(sqmahal(dist, m))
end

mahal_mee(dist::GaussVonMises, mat) = [mahalanobis(c, dist) for c in eachcol(mat)]

function save_mahal_mee!(d_q, name, i, model, mat)
    d_q[:, MEE_IDX[name], i] = quantile(mahal_mee(model, mat), QUANTILES)
end

function gvm_to_eucl_mc(gvm, N)
    # Package sampler returns a Vector of SVectors (one per sample)
    samples = reduce(hcat, rand(Random.default_rng(), gvm, (N,)))
    eucl = mapslices(mee_to_newton, samples, dims = 1)
    return fit(MvNormal, eucl)
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
    ref_samples = [SVector{L_STATE}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]
    current_ref_samples = copy(ref_samples)

    # Resampling to have a baseline
    alt_samples = [SVector{L_STATE}(col) for col in eachcol(rand(starting_dist, NSAMPLES_TRUTH))]
    current_alt_samples = copy(alt_samples)

    μ₀ = mean(starting_dist)
    P₀ = cov(starting_dist)

    μ₀_kepler = mean(starting_dist_kep)
    P₀_kepler = cov(starting_dist_kep)

    # Initial GVM (curvature Γ = 0) shared by all propagated GVM variants
    gvm0 = moment_match(starting_dist_kep, L_COORD)

    ts = Δt:Δt:t1
    nsteps = length(ts)

    h5open(outfile, "w") do f

        attrs(f)["NSAMPLES_TRUTH"] = NSAMPLES_TRUTH
        attrs(f)["NSAMPLES_MC_MAPPING"] = NSAMPLES_MC_MAPPING
        attrs(f)["delta_t"] = Δt
        attrs(f)["methods"] = map(v->v.name, METHODS)
        attrs(f)["mee_methods"] = MEE_METHODS
        attrs(f)["QUANTILES"] = QUANTILES
        f["ts"] = collect(ts)

        d_truth = create_dataset(
            f,
            "truth",
            Float32,
            (L_STATE, NSAMPLES_TRUTH, nsteps);
            chunk = (L_STATE, NSAMPLES_TRUTH, 1),
            deflate = 3,
        )

        d_means = create_dataset(
            f,
            "means",
            Float32,
            (L_STATE, length(METHODS), nsteps);
            chunk = (L_STATE, length(METHODS), 1),
            deflate = 3,
        )

        d_covsm = create_dataset(
            f,
            "covs",
            Float32,
            (L_STATE, L_STATE, length(METHODS), nsteps);
            chunk = (L_STATE, L_STATE, length(METHODS), 1),
            deflate = 3,
        )

        d_mahal_q = create_dataset(
            f,
            "mahalanobis_q",
            Float32,
            (length(QUANTILES), length(METHODS), nsteps),
        )

        d_mahal_mee = create_dataset(
            f,
            "mahalanobis_q_mee",
            Float32,
            (length(QUANTILES), length(MEE_METHODS), nsteps),
        )

        for i = 1:nsteps
            t = Δt * i
            current_ref_samples = run_monte_carlo(fm, current_ref_samples, Δt)
            current_alt_samples = run_monte_carlo(fm, current_alt_samples, Δt)
            current_ref_samples_kepler = map(newton_to_mee, current_ref_samples)

            mc_mat = reduce(hcat, current_ref_samples)
            mc_alt_mat = reduce(hcat, current_alt_samples)
            mc_mat_kepler = reduce(hcat, current_ref_samples_kepler)

            d_truth[:, :, i] = mc_mat

            # Slots that have no Cartesian-space meaning start as NaN
            d_mahal_q[:, :, i] = fill(Float32(NaN), length(QUANTILES), length(METHODS))
            d_mahal_mee[:, :, i] =
                fill(Float32(NaN), length(QUANTILES), length(MEE_METHODS))

            ut_dist = run_ut(fm, μ₀, P₀, t; equilibrate = true)
            stm_dist = run_stm(fm, μ₀, P₀, t; equilibrate = true)
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
                mee_to_newton,
                mean(ut_dist_kepler),
                cov(ut_dist_kepler),
                α = 1e-1,
                equilibrate = true,
            )
            stm_dist_kepler_eucl_ut = ut_propagate(
                mee_to_newton,
                mean(stm_dist_kepler),
                cov(stm_dist_kepler),
                α = 1e-1,
                equilibrate = true,
            )

            save_dist!(d_means, d_covsm, "UT (MEE->EUCL UT)", i, ut_dist_kepler_eucl_ut)
            save_dist!(d_means, d_covsm, "STM (MEE->EUCL UT)", i, stm_dist_kepler_eucl_ut)

            # MEE -> Euclid mapping via sampling
            samples_ut_dist_kepler = rand(ut_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_ut_dist_kepler_eucl =
                mapslices(mee_to_newton, samples_ut_dist_kepler, dims = 1)
            ut_dist_kepler_eucl_mc = fit(MvNormal, samples_ut_dist_kepler_eucl)

            save_dist!(d_means, d_covsm, "UT (MEE->EUCL MC)", i, ut_dist_kepler_eucl_mc)

            samples_stm_dist_kepler = rand(stm_dist_kepler, NSAMPLES_MC_MAPPING)
            samples_stm_dist_kepler_eucl =
                mapslices(mee_to_newton, samples_stm_dist_kepler, dims = 1)
            stm_dist_kepler_eucl_mc = fit(MvNormal, samples_stm_dist_kepler_eucl)
            save_dist!(d_means, d_covsm, "STM (MEE->EUCL MC)", i, stm_dist_kepler_eucl_mc)

            # Reference and alternate sampling MC
            ref_mc = fit(MvNormal, mc_mat)
            alt_mc = fit(MvNormal, mc_alt_mat)

            save_dist!(d_means, d_covsm, "Ref MC", i, ref_mc)
            save_dist!(d_means, d_covsm, "Resampled MC", i, alt_mc)

            # Control moment-matched GVM of the UT / STM MEE Gaussian (Γ = 0)
            gvm_mm_ut = moment_match(ut_dist_kepler, L_COORD)
            gvm_mm_stm = moment_match(stm_dist_kepler, L_COORD)

            name_mm = "UT moment-match GVM (MEE->EUCL MC)"
            gvm_mm_ut_eucl = gvm_to_eucl_mc(gvm_mm_ut, NSAMPLES_MC_MAPPING)
            save_dist!(d_means, d_covsm, name_mm, i, gvm_mm_ut_eucl)

            # Propagated GVMs
            mee_models = Pair{String,Any}[
                "UT (MEE)"=>ut_dist_kepler,
                "STM (MEE)"=>stm_dist_kepler,
                "UT (MEE) GVM moment-match"=>gvm_mm_ut,
                "STM (MEE) GVM moment-match"=>gvm_mm_stm,
            ]
            eucl_models = Pair{String,Any}[
                "UT"=>ut_dist,
                "STM"=>stm_dist,
                "UT (MEE->EUCL UT)"=>ut_dist_kepler_eucl_ut,
                "STM (MEE->EUCL UT)"=>stm_dist_kepler_eucl_ut,
                "UT (MEE->EUCL MC)"=>ut_dist_kepler_eucl_mc,
                "STM (MEE->EUCL MC)"=>stm_dist_kepler_eucl_mc,
                "Ref MC"=>ref_mc,
                "Resampled MC"=>alt_mc,
                name_mm=>gvm_mm_ut_eucl,
            ]

            for v in GVM_VARIANTS
                try
                    gvm_end =
                        run_gvm(fm_kepler, gvm0, t; slopes = v.slopes, refine = v.refine)
                    gvm_eucl = gvm_to_eucl_mc(gvm_end, NSAMPLES_MC_MAPPING)
                    save_dist!(d_means, d_covsm, v.name, i, gvm_eucl)
                    push!(mee_models, v.mee_name => gvm_end)
                    push!(eucl_models, v.name => gvm_eucl)
                catch e
                    @warn "GVM variant failed" variant = v.name t exception =
                        (e, catch_backtrace())
                    save_nan!(d_means, d_covsm, d_mahal_q, v.name, i)
                end
            end

            # Mahalanobis quantiles, Cartesian space (all methods vs. MC samples)
            for (name, model) in eucl_models
                save_mahal!(d_mahal_q, name, i, model, mc_mat)
            end

            # Mahalanobis quantiles, MEE space (MC samples converted to MEE)
            for (name, model) in mee_models
                save_mahal_mee!(d_mahal_mee, name, i, model, mc_mat_kepler)
            end

            @info "Propagation complete" t
            flush(f)
        end

    end

    return outfile
end

function main()
    END_T = 3600.0 * 24.0 * 4.0
    DELTA_T = 3600.0

    forces = DRAG_LOGβ_MEAN === nothing ? EARTH_FM_WITH_J2_TUPLE :
             (EARTH_FM_WITH_J2_TUPLE..., SimpleDragForce(DRAG_ρ₀, DRAG_r₀, DRAG_H))
    fm = ForceModel(GM_EARTH, forces, Val(true))
    fmk = ForceModel(GM_EARTH, forces, Val(false))

    μ = [6771e3, 0, 0, 0, 10142, 0]
    # μ = [9000e3, 0, 0, 0, 8620, 0]
    σ = Diagonal([1e3, 1e3, 1e3, 0.02, 0.02, 0.02])

    starting_dist = with_drag(MvNormal(μ, σ^2))

    # Transform to MEE coordinates, assuming normal after the non-linear transform
    starting_dist_kep = with_drag(ut_propagate(
        v -> euclid_to_mee(v..., GM_EARTH),
        μ,
        σ^2,
        α = 1e-1,
        equilibrate = true,
    ); mee = true)

    run_comparison(fm, fmk, starting_dist, starting_dist_kep, END_T, DELTA_T)
end

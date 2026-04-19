# PropagateGaussianOnMee.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox

include("Utils.jl")

fm_eucl = EARTH_FM_WITH_J2_NEWTON
fm_kepl = EARTH_FM_WITH_J2_KEPLER

μ = [9000e3, 0, 0, 0, 6620, 0]
σ = Diagonal([100e3, 10e3, 10e3, 1.0, 100.0, 1.0])

starting_dist = MvNormal(μ, σ^2)
starting_dist_mee = ut_propagate(
    v -> euclid_to_mee(v..., GM_EARTH),
    μ,
    σ^2,
)

function run_animation(fm_eucl, fm_kepl, starting_dist, starting_dist_mee, t1, Δt, filename="animation.mp4")
    starting_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, 1000))]
    starting_samples_mee = [SVector{6}(col) for col in eachcol(rand(starting_dist_mee, 1000))]
    
    current_samples = copy(starting_samples)
    current_samples_mee = copy(starting_samples_mee)
    μ₀ = mean(starting_dist_mee)
    P₀ = cov(starting_dist_mee)

    fig = Figure(size=(1366, 768))
    ax = Axis(fig[1,1])

    record(fig, filename, 0:Δt:t1, framerate=15, px_per_unit=2) do t
        t == 0 && return
        empty!(ax)

        current_samples = run_monte_carlo(fm_eucl, current_samples, Δt)
        current_samples_mee = run_monte_carlo(fm_kepl, current_samples_mee, Δt)
        # Ut always in MEE
        ut_dist = run_ut(fm, μ₀, P₀, t)


        mc_euclidean = stack(current_samples)
        scatter!(ax, mc_euclidean[1,:], mc_euclidean[2,:], color=(:red, 0.1), label="MC (MEE) 1000 samples")

        mee_mc = stack(current_samples_mee)
        mee_mc_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), mee_mc, dims=1)
        scatter!(ax, mee_mc_euclidean[1,:], mee_mc_euclidean[2,:], color=(:orange, 0.1), label="MC (Euclidean) 1000 samples")

        ut = rand(ut_dist, 1000)
        ut_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), ut, dims=1)
        scatter!(ax, ut_euclidean[1,:], ut_euclidean[2,:], color=(:blue, 0.1), label = "UT 1000 samples")
        

        r = 6000e3
        center = mean(mc_euclidean, dims=2)
        # center = [2000e3, 0e3, 0]
        xlims!(ax, center[1] - r, center[1] + r)
        ylims!(ax, center[2] - r, center[2] + r)

        axislegend(ax)
        @info t
    end
end

run_animation(fm_eucl, fm_kepl, starting_dist, starting_dist_mee, 7200, 20)

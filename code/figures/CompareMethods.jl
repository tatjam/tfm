# CompareMethods.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox

include("Utils.jl")

fm = EARTH_FM_WITH_J2_NEWTON

μ = [9000e3, 0, 0, 0, 6620, 0]
σ = Diagonal([100e3, 10e3, 10e3, 1.0, 100.0, 1.0])

starting_dist = MvNormal(μ, σ^2)

function run_animation(fm, starting_dist, t1, Δt, filename="animation.mp4")
    starting_samples = [SVector{6}(col) for col in eachcol(rand(starting_dist, 10000))]
    
    current_samples = copy(starting_samples)
    μ₀ = mean(starting_dist)
    P₀ = cov(starting_dist)

    fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1,1])

    record(fig, filename, 0:Δt:t1, framerate=15, px_per_unit=2) do t
        t == 0 && return
        empty!(ax)

        current_samples = run_monte_carlo(fm, current_samples, Δt)
        ut_dist = run_ut(fm, μ₀, P₀, t)
        stm_dist = run_stm(fm, μ₀, P₀, t)

        mc = stack(current_samples)
        scatter!(ax, mc[1,:], mc[2,:], color=(:orange, 0.05), label="MC")
        plot_ellipse(ut_dist, 1, 2, color=(:blue, 0.5), label="UT")
        plot_ellipse(stm_dist, 1, 2, color=(:red, 0.5), label="STM")

        # ut_samples = stack(rand(ut_dist, 10000))
        # stm_samples = stack(rand(stm_dist, 10000))
        # scatter!(ax, ut_samples[1,:], ut_samples[2,:], color=(:blue, 0.1))
        # scatter!(ax, stm_samples[1,:], stm_samples[2,:], color=(:red, 0.1))

        r = 10000e3
        center = mean(ut_dist)
        # center = [2000e3, 0e3, 0]
        xlims!(ax, center[1] - r, center[1] + r)
        ylims!(ax, center[2] - r, center[2] + r)

        axislegend(ax)
        @info t
    end
end

run_animation(fm, starting_dist, 7200, 20)

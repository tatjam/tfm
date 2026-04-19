# NaiveGaussianOnMee.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox

include("Utils.jl")

μ = [9000e3, 0, 0, 0, 6620, 0]
σ = Diagonal([2000e3, 10e3, 10e3, 1.0, 100.0, 1.0])

starting_dist = MvNormal(μ, σ^2)

function mee_and_back(dist)

    starting_dist_mee = ut_propagate(
        v -> euclid_to_mee(v..., GM_EARTH),
        μ,
        σ^2,
        α=1e-1
    )

    starting_dist_back = ut_propagate(
        v -> mee_to_euclid(v...,GM_EARTH),
        mean(starting_dist_mee),
        cov(starting_dist_mee)
    )

    fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1,1])


    i = 1
    j = 2
    # plot_ellipse(starting_dist, i, j, color=(:red, 1))

    samples_mee = rand(starting_dist_mee, 10000)
    # scatter!(ax, samples_mee[i,:], samples_mee[j,:], color=(:orange, 0.05))
    results = SVector{6,Float64}[]
    for col in eachcol(samples_mee)
        try
            push!(results, mee_to_euclid(col..., GM_EARTH))
        catch
        end
    end
    samples_mee_back = reduce(hcat, results)

    scatter!(ax, samples_mee_back[i,:], samples_mee_back[j,:], color=(:orange, 0.05))

    plot_ellipse(starting_dist_back, i, j, color=(:blue, 1))

    fig
end

mee_and_back(starting_dist)


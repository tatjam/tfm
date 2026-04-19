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

fm = EARTH_FM_WITH_J2_KEPLER

μ = [9000e3, 0, 0, 0, 6620, 0]
σ = Diagonal([1e3, 1e3, 1e3, 1.0, 1.0, 1.0])

starting_dist = MvNormal(μ, σ^2)

starting_dist_mee = ut_propagate(
    v -> kepler_to_mee(euclid_to_kepler(v..., GM_EARTH)...),
    μ,
    σ^2
)

starting_dist_back = ut_propagate(
    v -> kepler_to_euclid(mee_to_kepler(v...)..., GM_EARTH),
    mean(starting_dist_mee),
    nearest_pd_matrix(cov(starting_dist_mee))
)

fig = Figure(size=(1920, 1080))
ax = Axis(fig[1,1])

plot_ellipse(starting_dist, 1, 2, color=(:red, 0.5))
plot_ellipse(starting_dist_back, 1, 2, color=(:blue, 0.5))

fig


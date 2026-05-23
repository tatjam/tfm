# PropagateGVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------


using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox

include("Utils.jl")

fm_kepl = EARTH_FM_WITH_J2_KEPLER

orbit_u0_mee = [7136.635, 0, 0, 0, 0] * 1e3
P = Diagonal([(20e3)^2, 10e-6, 10e-6, 10e-6, 10e-6])
β = [0.0, 0.0, 0.0, 0.0, 0.0]
Γ = zeros(5, 5)
κ = 3.282806e7

dist = GaussVonMises(orbit_u0_mee, 0.0, β, Γ, κ, P=P)

function run_animation(fm, starting_dist, t1, Δt, filename="animation.mp4")
    starting_samples = [SVector{6}(col) for col in eachcol(reduce(hcat, rand(starting_dist, 10000)))]
    
    current_samples = copy(starting_samples)

    fig = Figure(size=(1366, 768))
    ax = Axis(fig[1,1])
    ax_state = Axis(fig[1,2])
    state_x = 1
    state_y = 6

    record(fig, filename, 0:Δt:t1, framerate=15, px_per_unit=2) do t
        t == 0 && return
        empty!(ax)
        empty!(ax_state)

        current_samples = run_monte_carlo(fm, current_samples, Δt)
        gvm_dist = run_gvm(fm, starting_dist, t)

        mc = stack(current_samples)
        scatter!(ax_state, mc[state_x,:], mc[state_y,:], color=(:orange, 0.1))
        mc_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), mc, dims=1)
        scatter!(ax, mc_euclidean[1,:], mc_euclidean[2,:], color=(:orange, 0.1), label="MC 1000 samples")

        gvm = reduce(hcat, rand(gvm_dist, 10000))
        scatter!(ax_state, gvm[state_x,:], gvm[state_y,:], color=(:blue, 0.1))
        gvm_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), gvm, dims=1)
        scatter!(ax, gvm_euclidean[1,:], gvm_euclidean[2,:], color=(:blue, 0.1), label = "GVM 1000 samples")
        

        r = 6000e3
        center = mean(mc_euclidean, dims=2)
        # center = [2000e3, 0e3, 0]
        xlims!(ax, center[1] - r, center[1] + r)
        ylims!(ax, center[2] - r, center[2] + r)

        axislegend(ax)
        @info t
    end
end

run_animation(fm_kepl, dist, 7000, 50)


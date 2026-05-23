# PropagateGVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------


using OrbitalUncertainty
using GLMakie
using Distributions
using StaticArrays
using LinearAlgebra
using SatelliteToolbox
using Printf

include("Utils.jl")

fm_kepl = EARTH_FM_WITH_J2_KEPLER

orbit_u0_mee = [7136.635, 0, 0, 0, 0] * 1e3
P = Diagonal([(20e3)^2, 10e-6, 10e-6, 10e-6, 10e-6])
# P = Diagonal([20e3, 10e-6, 10e-6, 10e-6, 10e-6])
β = [0.0, 0.0, 0.0, 0.0, 0.0]
Γ = zeros(5, 5)
κ = 3.282806e7

dist = GaussVonMises(orbit_u0_mee, 0.0, β, Γ, κ, P=P)

function run_animation(fm, starting_dist, t1, Δt, filename="animation.mp4")
    starting_samples = [SVector{6}(col) for col in eachcol(reduce(hcat, rand(starting_dist, 10000)))]
    
    current_samples = copy(starting_samples)
    fig = Figure(size=(1366*2, 768))
    ax = Axis(fig[1,1], title="Euclidean space")
    ax_state = Axis(fig[1,2], title="MEE space")
    ax_beta = Axis(fig[1,3], title="β")
    ax_gamma = Axis(fig[1,4], title="Γ")
    state_x = 1
    state_y = 6
    record(fig, filename, 0:Δt:t1, framerate=15, px_per_unit=2) do t
        t == 0 && return
        empty!(ax)
        empty!(ax_state)
        empty!(ax_beta)
        empty!(ax_gamma)
        current_samples = run_monte_carlo(fm, current_samples, Δt)
        gvm_dist = run_gvm(fm, starting_dist, t)
        mc = stack(current_samples)
        scatter!(ax_state, mc[state_x,:], mc[state_y,:], color=(:orange, 0.1))
        mc_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), mc, dims=1)
        scatter!(ax, mc_euclidean[1,:], mc_euclidean[2,:], color=(:orange, 0.1), label="MC 1000 samples")
        gvm = reduce(hcat, rand(gvm_dist, 10000))
        scatter!(ax_state, gvm[state_x,:], gvm[state_y,:], color=(:blue, 0.1))
        gvm_euclidean = mapslices(v -> mee_to_euclid.(v..., GM_EARTH), gvm, dims=1)
        scatter!(ax, gvm_euclidean[1,:], gvm_euclidean[2,:], color=(:blue, 0.1), label="GVM 1000 samples")
        
        # β heatmap - it's a vector so show as a column
        β = gvm_dist.β
        heatmap!(ax_beta, reshape(β, :, 1), colormap=:RdBu, colorrange=(-maximum(abs.(β)) - 1e-10, maximum(abs.(β)) + 1e-10))
        text!(ax_beta, [Point2f(0, i-1) for i in 1:length(β)], 
              text=[@sprintf("%.2e", β[i]) for i in 1:length(β)],
              align=(:center, :center), fontsize=10)

        # Γ heatmap - it's an n×n matrix
        Γ = gvm_dist.Γ
        heatmap!(ax_gamma, Γ, colormap=:RdBu, colorrange=(-maximum(abs.(Γ)) - 1e-10, maximum(abs.(Γ)) + 1e-10))
        for i in 1:size(Γ, 1), j in 1:size(Γ, 2)
            text!(ax_gamma, Point2f(i, j), 
                  text=@sprintf("%.2e", Γ[i,j]),
                  align=(:center, :center), fontsize=8)
        end

        r = 6000e3
        center = mean(mc_euclidean, dims=2)
        xlims!(ax, center[1] - r, center[1] + r)
        ylims!(ax, center[2] - r, center[2] + r)
        axislegend(ax)
        @info "t=$t β_max=$(maximum(abs.(β))) Γ_max=$(maximum(abs.(Γ)))"
    end
end
run_animation(fm_kepl, dist, 48000, 600)


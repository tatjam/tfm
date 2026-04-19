# Utils.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using GLMakie

function plot_ellipse(dist::MvNormal, i, j; nσ=2, n=200, color=(:blue, 1.0), label="")
    # Extract the marginal statistics
    μ = reshape(mean(dist)[[i, j]], 2, 1)
    P = cov(dist)[[i, j], [i, j]]

    λ, v = eigen(P)

    # Build the ellipse
    θ = range(0, 2π, n)
    circle = [cos.(θ), sin.(θ)]
    ellipse = stack(v * Diagonal(nσ .* sqrt.(λ)) * circle)' .+ μ

    lines!(ellipse[1,:], ellipse[2,:], color=color, label=label)
end

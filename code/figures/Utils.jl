# Utils.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using GLMakie
using Distributions

function plot_ellipse(
    dist::MvNormal,
    i,
    j;
    nσ = 2,
    n = 200,
    color = (:blue, 1.0),
    label = "",
)
    # Extract the marginal statistics
    μ = reshape(mean(dist)[[i, j]], 2, 1)
    P = cov(dist)[[i, j], [i, j]]

    λ, v = eigen(P)

    # Build the ellipse
    θ = range(0, 2π, n)
    circle = [cos.(θ), sin.(θ)]
    ellipse = stack(v * Diagonal(nσ .* sqrt.(λ)) * circle)' .+ μ

    lines!(ellipse[1, :], ellipse[2, :], color = color, label = label)
end

function plot_earth(ax3)
    R_EARTH_KM = 6371.0
    mesh!(ax3, Sphere(Point3f(0, 0, 0), R_EARTH_KM), color = (:steelblue, 0.6))
end

function plot_rv_ellipsoid(
    ax3,
    μ,
    Σ;
    nσ = 3,
    resolution = 30,
    is = 1:3,
    scale = 1e-3,
    relative = false,
    R = nothing,
)
    r = μ[is] .* scale
    Σ_rr = Σ[is, is] .* (scale ^ 2)

    c = relative ? zeros(is) : r

    # Base sphere
    u = range(0, 2π, resolution)
    v = range(-π/2, π/2, resolution)

    # This builds a resolution x resolution grid of values
    xs = cos.(u)' .* cos.(v)
    ys = sin.(u)' .* cos.(v)
    zs = ones(length(u))' .* sin.(v)

    # We flatten the points to vector
    xs = vec(xs)
    ys = vec(ys)
    zs = vec(zs)

    # Finally, stack all the points, note the transpose to end up with
    # a 3xN matrix, i.e. a matrix of column vectors
    pts = hcat(xs, ys, zs)'

    # Rotate the xyz axes to the eigen vectors, scaled by their roots and nσ
    λ, λv = eigen(Σ_rr)
    λ = max.(λ, 0)
    tform = λv * Diagonal(nσ .* sqrt.(λ))

    pts = (tform * pts) .+ c

    # Bring back the grid of points for surface
    xs = reshape(pts[1, :], resolution, resolution)
    ys = reshape(pts[2, :], resolution, resolution)
    zs = reshape(pts[3, :], resolution, resolution)

    surface!(
        ax3,
        xs,
        ys,
        zs;
        color = fill(1, size(xs)),
        colormap = [:steelblue],
        transparency = true,
    )
    wireframe!(ax3, xs, ys, zs; color = (:gray), linewidth = 1)

    if !isnothing(R)
        L = 1.5 * maximum(nσ .* sqrt.(λ))
        origin = Point3f(c)
        cols = [:red, :green, :blue]
        for (k, (lbl, col)) in enumerate(zip(["R", "T", "N"], cols))
            tip = Point3f(origin + L * Point3f(R[1:3, k]))
            lines!(ax3, [origin, tip]; color = col, linewidth = 1)
            text!(ax3, tip; text = lbl, color = col, fontsize = 16)
        end
    end

end

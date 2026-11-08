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

function plot_ellipsoid_and_displaced!(
    ax3,
    Σ,
    R;
    nσ = 3,
    resolution = 15,
    scale = 1e-3,
    dt = 600,
    colormap = :viridis,
    alpha = 0.5,
)
    Σ_rr = Σ[1:3, 1:3] .* (scale ^ 2)
    Σ_vr = Σ[4:6, 1:3] .* (scale ^ 2)
    Σ_vv = Σ[4:6, 4:6] .* (scale ^ 2)

    # Check if we have any position-scale correlation, to avoid rounding
    # error from triggering divide by near zero
    σr = sqrt.(diag(Σ_rr))
    σv = sqrt.(diag(Σ_vv))
    negligible = maximum(abs.(Σ_vr ./ (σv * σr'))) < 1e-6

    # position offset -> expected velocity offset [1/s]
    K = negligible ? zeros(3, 3) : Σ_vr / Σ_rr

    λ, λv = eigen(Symmetric(Σ_rr))
    λ = max.(λ, 0)
    radii = nσ .* sqrt.(λ)
    tform = λv * Diagonal(radii)

    # Base sphere
    u = range(0, 2π, resolution)
    v = range(-π/2, π/2, resolution)
    xs = vec(cos.(u)' .* cos.(v))
    ys = vec(sin.(u)' .* cos.(v))
    zs = vec(ones(length(u))' .* sin.(v))
    pts = hcat(xs, ys, zs)'

    # Starting covariance ellipsoid
    drs = tform * pts

    # Full position covariance after dt (without dynamics)
    Σ2 = Symmetric(Σ_rr + dt * (Σ_vr + Σ_vr') + dt^2 * Σ_vv)
    λ2, λv2 = eigen(Σ2)
    tform2 = λv2 * Diagonal(nσ .* sqrt.(max.(λ2, 0)))
    drs2 = tform2 * pts

    function draw(P, wirecolor, color)
        X = reshape(P[1, :], resolution, resolution)
        Y = reshape(P[2, :], resolution, resolution)
        Z = reshape(P[3, :], resolution, resolution)
        surface!(
            ax3,
            X,
            Y,
            Z;
            color = color,
            shading = NoShading,
            transparency = true,
            rasterize = 1,
        )
        wireframe!(ax3, X, Y, Z; color = (wirecolor, 0.5), linewidth = 0.7)
    end

    draw(drs, :black, :slateblue)
    draw(drs2, :red, :tomato2)

    # Legend
    t = dt >= 3600 ? "$(round(Int, dt / 3600)) h" : "$(round(Int, dt / 60)) min"
    axislegend(
        ax3,
        [LineElement(color = :black), LineElement(color = :red)],
        ["Original", "Tras $t"];
        position = :rt,
    )

    # Frame length (covers both ellipsoids)
    L = 1.3 * max(maximum(norm.(eachcol(drs))), maximum(norm.(eachcol(drs2))))

    # Cover both ellipsoids and the frame tips + labels
    ext = [
        max(
            maximum(abs.(drs[k, :])),
            maximum(abs.(drs2[k, :])),
            L * maximum(abs.(R[k, :])),
        ) for k = 1:3
    ]
    ext .*= 1.1
    lo = - ext
    hi = ext

    # Fix the limits so the walls are exactly where we draw the projections
    xlims!(ax3, lo[1], hi[1])
    ylims!(ax3, lo[2], hi[2])
    zlims!(ax3, lo[3], hi[3])

    # Walls: x and y at the upper limit, z at the lower limit
    # (This is the default Makie view)
    wall = [hi[1], hi[2], lo[3]]

    # RTN frame vectors (columns of R are the axes)
    origin = Point3f(zeros(3))
    for (k, (lbl, col)) in enumerate(zip(["R", "T", "N"], [:red, :green, :blue]))
        tip = Point3f(origin + L * Point3f(R[1:3, k]))
        lines!(ax3, Point3f[origin, tip]; color = col, linewidth = 1, overdraw = true)
        text!(ax3, tip; text = lbl, color = col, fontsize = 16, overdraw = true)
    end

    # Ellipsoid projections on the walls
    θ = range(0, 2π, 100)
    circ = [cos.(θ)'; sin.(θ)']
    for (M, col) in ((tform, :black), (tform2, :red))
        for k = 1:3
            keep = setdiff(1:3, k)
            M2 = M[keep, :]
            λp, Vp = eigen(Symmetric(M2 * M2'))
            e = Vp * Diagonal(sqrt.(max.(λp, 0))) * circ
            P = zeros(3, size(e, 2))
            P[keep, :] = e
            P[k, :] .= wall[k]
            lines!(
                ax3,
                [Point3f(P[:, i]) for i = 1:size(P, 2)];
                color = (col, 0.8),
                linewidth = 1.5,
            )
        end
    end

    # RTN axis projections, also on the walls
    for (j, col) in enumerate([:red, :green, :blue]), k = 1:3
        p0 = Vector{Float64}(zeros(3))
        p1 = p0 .+ L .* R[1:3, j]
        p0[k] = wall[k]
        p1[k] = wall[k]
        lines!(
            ax3,
            Point3f[Point3f(p0), Point3f(p1)];
            color = col,
            linewidth = 1,
            linestyle = :dash,
        )
    end

    return dt
end

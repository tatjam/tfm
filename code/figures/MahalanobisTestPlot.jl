# MahalanobisTestPlot.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using GLMakie
using HDF5
using Statistics

include("MahalanobisTest.jl")

METHOD_PLOT = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DIST_LOG = true

function plot_band!(ax, t, q, m)
    b  = band!(ax, t, sqrt.(q[1, :]), sqrt.(q[3, :]), color = (m.color, 0.25))
    md = lines!(ax, t, sqrt.(q[2, :]), color = m.color)
    lo = lines!(ax, t, sqrt.(q[1, :]), color = m.color, linestyle = :dash)
    hi = lines!(ax, t, sqrt.(q[3, :]), color = m.color, linestyle = :dash)
    return [b, md, lo, hi]
end

function main_plot()
    @info "Loading"
    file = get(ARGS, 1, "comparison.h5")
    ts, mahalanobis_q = h5open(file, "r") do f
        read(f["ts"]), read(f["mahalanobis_q"])
    end

    ts /= 3600.0


    plot_methods = [(i, m) for (i, m) in enumerate(METHODS) if i in METHOD_PLOT]

    GLMakie.activate!()
    fig = Figure(size = (1800, 900))
    @info "Mahalanobis Plot"

    @info "Mahalanobis quartiles Plot"
    ax_q = Axis(
        fig[1, 1],
        xlabel = "t [h]",
        ylabel = "Mahalanobis (quartiles)",
        yscale = DIST_LOG ? log10 : sqrt,
        yautolimitmargin = (0.0, 0.05)
    )


    panel = GridLayout(fig[1, 2], tellheight = false)
    
    for (row, (i, m)) in enumerate(plot_methods)
        band_plots = plot_band!(ax_q, ts, mahalanobis_q[:, i, :], m)

        toggle = Toggle(fig, active = true, buttoncolor = m.color)
        on(toggle.active; update = true) do v
            for p in band_plots
                p.visible = v
            end
        end
        panel[row, 1] = toggle
        panel[row, 2] = Label(fig, m.name; halign = :left)
    end

    display(fig)

    return fig

end

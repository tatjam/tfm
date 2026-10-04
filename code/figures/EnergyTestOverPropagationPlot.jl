# EnergyTestOverPropagationPlot.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using GLMakie
using HDF5
using Statistics

# Methods stored in the HDF5 file, in sample order ("mc,ut,stm,ut_kep,stm_kep").
# Toggle which ones to visualize by editing METHOD_PLOT.
METHODS = [
    (name = "MC", color = :darkorange),
    (name = "UT", color = :seagreen),
    (name = "STM", color = :purple),
    (name = "UT (MEE)", color = :royalblue),
    (name = "STM (MEE)", color = :crimson),
]
METHOD_PLOT = [1, 2, 3, 4, 5]

# Plot the energy curves on a logarithmic y-axis.
ENERGY_LOG = true

function main()
    file = get(ARGS, 1, "comparison.h5")
    ts, energies, truth, samples = h5open(file, "r") do f
        read(f["ts"]), read(f["energies"]), read(f["truth"]), read(f["samples"])
    end
    # truth:   (6, NSAMPLES_TRUTH, nsteps)  rows 1:3 = position [m]
    # samples: (6, NSAMPLES, 5, nsteps)     rows 1:3 = position [m]

    nsteps = length(ts)

    plot_methods = [(i, m) for (i, m) in enumerate(METHODS) if i in METHOD_PLOT]

    # Single window: energy plot + time slider (left), 3D orbital view (middle),
    # velocity vector scatter at the selected time (right)
    GLMakie.activate!()
    fig = Figure(size = (1800, 900))
    ax = Axis(
        fig[1, 1],
        xlabel = "t [s]",
        ylabel = "Energy",
        yscale = ENERGY_LOG ? log10 : identity,
    )
    for (i, m) in plot_methods
        lines!(ax, ts, energies[i, :], color = m.color, label = m.name)
    end
    axislegend(ax)

    # Slider picking the time visualized in the orbital view
    sg = SliderGrid(
        fig[2, 1],
        (label = "t", range = 1:nsteps, format = i -> "$(round(ts[i] / 86400, digits = 2)) d"),
    )
    i_time = sg.sliders[1].value
    set_close_to!(sg.sliders[1], nsteps)
    vlines!(ax, lift(i -> ts[i], i_time), color = :gray, linewidth = 1)

    # 3D orbital view at the selected time
    ax3 = Axis3(
        fig[1:2, 2],
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [km]",
        title = lift(
            i -> "Orbital view @ t = $(round(ts[i] / 86400, digits = 2)) days",
            i_time,
        ),
        aspect = :data,
        xgridcolor = :gray,
        ygridcolor = :gray,
        zgridcolor = :gray,
    )

    # Ground-truth Monte Carlo cloud (backdrop)
    scatter!(
        ax3,
        lift(i -> truth[1, :, i] ./ 1e3, i_time),
        lift(i -> truth[2, :, i] ./ 1e3, i_time),
        lift(i -> truth[3, :, i] ./ 1e3, i_time),
        color = (:gray, 0.3),
        label = "Truth MC",
    )

    # Same methods as the line plot
    for (i, m) in plot_methods
        scatter!(
            ax3,
            lift(j -> samples[1, :, i, j] ./ 1e3, i_time),
            lift(j -> samples[2, :, i, j] ./ 1e3, i_time),
            lift(j -> samples[3, :, i, j] ./ 1e3, i_time),
            color = (m.color, 0.6),
            label = m.name,
        )
    end

    # Mean trajectory of the truth cloud over the whole propagation
    mean_traj = dropdims(mean(truth, dims = 2), dims = 2)   # (6, nsteps)
    lines!(
        ax3,
        mean_traj[1, :] ./ 1e3,
        mean_traj[2, :] ./ 1e3,
        mean_traj[3, :] ./ 1e3,
        color = :gray,
        linewidth = 2,
        label = "Mean trajectory",
    )

    # Earth for context
    R_EARTH_KM = 6371.0
    mesh!(ax3, Sphere(Point3f(0, 0, 0), R_EARTH_KM), color = (:steelblue, 0.6))

    # Axis limits: 10% beyond the mean trajectory or Earth radius, whichever
    # is larger per axis, so outliers don't zoom everything out
    for (setlims, row) in zip((xlims!, ylims!, zlims!), 1:3)
        lim = max(1.1 * maximum(abs, mean_traj[row, :]) / 1e3, R_EARTH_KM)
        setlims(ax3, -lim, lim)
    end

    axislegend(ax3)

    # Velocity vector scatter at the selected time
    axv = Axis3(
        fig[1:2, 3],
        xlabel = "vx [km/s]",
        ylabel = "vy [km/s]",
        zlabel = "vz [km/s]",
        title = lift(
            i -> "Velocity vectors @ t = $(round(ts[i] / 86400, digits = 2)) days",
            i_time,
        ),
        aspect = :data,
        xgridcolor = :gray,
        ygridcolor = :gray,
        zgridcolor = :gray,
    )

    # Ground-truth Monte Carlo cloud (backdrop)
    scatter!(
        axv,
        lift(i -> truth[4, :, i] ./ 1e3, i_time),
        lift(i -> truth[5, :, i] ./ 1e3, i_time),
        lift(i -> truth[6, :, i] ./ 1e3, i_time),
        color = (:gray, 0.3),
        label = "Truth MC",
    )

    # Same methods as the line plot
    for (i, m) in plot_methods
        scatter!(
            axv,
            lift(j -> samples[4, :, i, j] ./ 1e3, i_time),
            lift(j -> samples[5, :, i, j] ./ 1e3, i_time),
            lift(j -> samples[6, :, i, j] ./ 1e3, i_time),
            color = (m.color, 0.6),
            label = m.name,
        )
    end

    # Axis limits: 10% beyond the 99.5% quantile of |component| over all
    # clouds, so outliers don't zoom everything out
    vlim = 1.1 * quantile(
        abs.(vcat(vec(truth[4:6, :, :]), vec(samples[4:6, :, :, :]))),
        0.995,
    ) / 1e3
    for setlims in (xlims!, ylims!, zlims!)
        setlims(axv, -vlim, vlim)
    end

    axislegend(axv)

    display(fig)

    return fig
end

main()

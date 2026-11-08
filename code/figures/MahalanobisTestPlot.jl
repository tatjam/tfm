# MahalanobisTestPlot.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using GLMakie
using HDF5
using Statistics

include("MahalanobisTest.jl")

# :eucl  only the Cartesian-space figure
# :mee   only the MEE-space figure
# :both  both, two separate windows
PLOT_MODE  = :mee
DIST_LOG   = false

# Which methods to include in each figure (indices into METHODS / MEE_METHODS)
const EUCL_METHODS = 1:length(METHODS)
const MEE_METHODS_PLOT = 1:length(MEE_METHODS)

const MEE_METHOD_COLORS = Dict(
    "UT (MEE)"                   => :deepskyblue,
    "STM (MEE)"                  => :crimson,
    "UT (MEE) GVM moment-match"  => :goldenrod,
    "STM (MEE) GVM moment-match" => :purple,
    "GVM :ad"                    => :forestgreen,
    "GVM :nodes"                 => :olivedrab,
    "GVM :ad refine"             => :darkgreen,
    "GVM :nodes refine"          => :limegreen,
)

function plot_band!(ax, t, q, color, label)
    q05 = sqrt.(q[1, :])
    q50 = sqrt.(q[2, :])
    q95 = sqrt.(q[3, :])
    b  = band!(ax, t, q05, q95, color = (color, 0.25))
    md = lines!(ax, t, q50, color = color, label = label)
    lo = lines!(ax, t, q05, color = color, linestyle = :dash)
    hi = lines!(ax, t, q95, color = color, linestyle = :dash)
    return [b, md, lo, hi]
end

function make_figure(; title, methods, data, ts)
    fig = Figure(size = (1500, 1100))   # a bit taller helps too

    ax = Axis(
        fig[1, 1],
        xlabel = "t [h]",
        ylabel = "√Mahalanobis (quantiles)",
        yscale = DIST_LOG ? log10 : identity,
        yautolimitmargin = (0.0, 0.05),
        title = title,
        aspect = nothing,
    )

    panel = GridLayout(
        fig[1, 2];
        tellheight = false,
        tellwidth  = true,
        valign     = :top,
    )
    colsize!(fig.layout, 2, Fixed(280))

    # Row 1 takes the full figure height
    rowsize!(fig.layout, 1, Relative(1.0))

    for (row, (i, name, color)) in enumerate(methods)
        ps = plot_band!(ax, ts, data[:, i, :], color, name)
        toggle = Toggle(fig, active = true, buttoncolor = color)
        on(toggle.active; update = true) do v
            for p in ps
                p.visible = v
            end
        end
        panel[row, 1] = toggle
        panel[row, 2] = Label(fig, name; halign = :left)
    end

    return fig
end

function make_eucl_plot(ts, mahal_q)
    methods = [(i, METHODS[i].name, METHODS[i].color) for i in EUCL_METHODS]
    return make_figure(;
        title   = "Cartesian space",
        methods = methods,
        data    = mahal_q,
        ts      = ts,
    )
end

function make_mee_plot(ts, mahal_q_mee)
    methods = [
        (i, MEE_METHODS[i], get(MEE_METHOD_COLORS, MEE_METHODS[i], :black))
        for i in MEE_METHODS_PLOT
    ]
    return make_figure(;
        title   = "MEE space",
        methods = methods,
        data    = mahal_q_mee,
        ts      = ts,
    )
end

function main_plot(; file = get(ARGS, 1, "comparison.h5"))
    @info "Loading" file
    ts, mahal_q, mahal_q_mee = h5open(file, "r") do f
        t  = read(f["ts"])
        qe = read(f["mahalanobis_q"])
        qm = haskey(f, "mahalanobis_q_mee") ? read(f["mahalanobis_q_mee"]) : nothing
        (t, qe, qm)
    end
    ts ./= 3600.0

    if !isnothing(mahal_q_mee) == false
        @warn "No `mahalanobis_q_mee` in file" file
    end

    GLMakie.activate!()

    figs = Dict{Symbol,Figure}()

    if PLOT_MODE in (:eucl, :both)
        figs[:eucl] = make_eucl_plot(ts, mahal_q)
        display(figs[:eucl])
    end

    if PLOT_MODE in (:mee, :both)
        if isnothing(mahal_q_mee)
            @warn "PLOT_MODE=$PLOT_MODE requested MEE plot but data is missing"
        else
            figs[:mee] = make_mee_plot(ts, mahal_q_mee)
            display(figs[:mee])
        end
    end

    return figs
end

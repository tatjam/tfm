# EnergyTestOverPropagationPlot.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using CairoMakie
using JLD2


function main()
    CairoMakie.activate!(px_per_unit=1)
    Makie.inline!(true)
    set_theme!(theme_dark())
    
    file = get(ARGS, 1, "comparison.jld2")
    ts, truth, sampled_states, energies = load(file, "ts", "truth", "sampled_states", "energies")

    fig = Figure(size=(1600, 900))
    ax = Axis(fig[1, 1])
    energy_names = ["MC", "UT", "STM", "UT (MEE)", "STM (MEE)"]
    for i in [1, 4, 5]
        lines!(ax, ts, energies[i,:], label = energy_names[i])
    end

    axislegend(ax)

    return fig
    
end

main()

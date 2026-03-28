# OrbitalUncertainty.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
module OrbitalUncertainty

using SatelliteToolbox
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using Distributions

# ForceModel
include("dynamics/ForceModel.jl")
export ForceModel, TwoBodyForce, J2Force
export acceleration, newton_model, propagate_orbit

include("propagators/MonteCarlo.jl")
export run_monte_carlo

include("propagators/UT.jl")
export SigmaVectors, run_ut

include("propagators/STM.jl")

# Utils 
EARTH_FM_WITH_J2 = ForceModel((TwoBodyForce(GM_EARTH), J2Force(GM_EARTH, EARTH_EQUATORIAL_RADIUS, EGM_1996_J2)))
export EARTH_FM_WITH_J2

end # module OrbitalUncertainty

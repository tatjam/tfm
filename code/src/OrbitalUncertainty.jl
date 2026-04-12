# OrbitalUncertainty.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
module OrbitalUncertainty

using SatelliteToolbox
using ReferenceFrameRotations
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using Distributions

# ForceModel
include("dynamics/ForceModel.jl")
include("dynamics/utils.jl")
export ForceModel, TwoBodyForce, J2Force
export acceleration, force_model, propagate_orbit
export kepler_to_mee, mee_to_kepler, kepler_to_array, kepler_to_euclid, euclid_to_kepler, isapprox_angle

include("propagators/MonteCarlo.jl")
export run_monte_carlo

include("propagators/UT.jl")
export SigmaVectors, ut_propagate, run_ut

include("propagators/STM.jl")
export run_stm

# Utils 
EARTH_FM_TUPLE = (TwoBodyForce(GM_EARTH),)
EARTH_FM_NEWTON = ForceModel(EARTH_FM_TUPLE, Val(true))
EARTH_FM_KEPLER = ForceModel(EARTH_FM_TUPLE, Val(false))
export EARTH_FM_NEWTON, EARTH_FM_KEPLER

EARTH_FM_WITH_J2_TUPLE = (TwoBodyForce(GM_EARTH), J2Force(GM_EARTH, EARTH_EQUATORIAL_RADIUS, EGM_1996_J2))
EARTH_FM_WITH_J2_NEWTON = ForceModel(EARTH_FM_WITH_J2_TUPLE, Val(true))
EARTH_FM_WITH_J2_KEPLER = ForceModel(EARTH_FM_WITH_J2_TUPLE, Val(false))
export EARTH_FM_WITH_J2_NEWTON, EARTH_FM_WITH_J2_KEPLER

end # module OrbitalUncertainty

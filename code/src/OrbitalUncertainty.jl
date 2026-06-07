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
using Bessels
using Random
using NLLSsolver
using Static: static

# ForceModel
include("dynamics/ForceModel.jl")
include("dynamics/Utils.jl")
export ForceModel, J2Force, EGM96Force
export acceleration, force_model, propagate_orbit
export kepler_to_euclid, kepler_to_array, euclid_to_kepler
export kepler_to_mee, mee_to_kepler
export euclid_to_mee, mee_to_euclid
export get_csn_basis, csn_acceleration_to_mee
export isapprox_angle

include("propagators/MonteCarlo.jl")
export run_monte_carlo

include("propagators/UT.jl")
export SigmaVectors, ut_propagate, run_ut, nearest_pd_matrix

include("propagators/STM.jl")
export run_stm

include("statistics/GVM.jl")
export GaussVonMises, decanonicalize, mahalanobis, canon_mahalanobis
include("propagators/GVM.jl")
export GVMSigmaVectors, b12, gvm_propagate, run_gvm

# Utils 
EARTH_FM_NEWTON = ForceModel(GM_EARTH, (), Val(true))
EARTH_FM_KEPLER = ForceModel(GM_EARTH, (), Val(false))
export EARTH_FM_NEWTON, EARTH_FM_KEPLER

EARTH_FM_WITH_J2_TUPLE = (J2Force(EARTH_EQUATORIAL_RADIUS, EGM_1996_J2),)
EARTH_FM_WITH_J2_NEWTON = ForceModel(GM_EARTH, EARTH_FM_WITH_J2_TUPLE, Val(true))
EARTH_FM_WITH_J2_KEPLER = ForceModel(GM_EARTH, EARTH_FM_WITH_J2_TUPLE, Val(false))
export EARTH_FM_WITH_J2_NEWTON, EARTH_FM_WITH_J2_KEPLER

end # module OrbitalUncertainty

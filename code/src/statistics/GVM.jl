# GVM.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Gauss Von-Mises distribution implementation, a natural representation
# for uncertainty in MEE coordinates. We use notation as defined in
# "Gauss von Mises Distribution for Improved Uncertainty Realism in
#  Space Situational Awareness", Joshua T. Horwood and Aubrey B. Poore, 2014.
# 
# For pure keplerian coordinates (a, e, i, ω, Ω, ν), this distribution is also
# appropriate, but careful choice of i = ω = Ω = 0 has to be done if these
# coordinates are to be considered Gaussian.

# STMEuclideanToMeeMap.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
include("./Cases.jl")

function main()
    for case_i = 1:length(CASES_MEAN_MEE)
        μ = CASES_MEAN_MEE[case_i]
        Σ = case_covariance(μ, CASES_COV_RTN_SCALED[case_i])

        @info μ, Σ
    end
end


main()

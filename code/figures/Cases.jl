# Cases.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using OrbitalUncertainty
using LinearAlgebra
using SatelliteToolbox
using Distributions

include("Utils.jl")

CASES_MEAN_MEE = [
    [6878e3, 0.001, 0.0, 0.471, 0.0, 0],
    [12028e3, 0, -0.74, 0.635, 0.0, 4.712],
    [12028e3, 0, -0.74, 0.635, 0.0, 1.57],
]

# Cases in RTN coordinates, scaled by the radius and velocity of starting point
#! format: off
CASES_COV_RTN_SCALED = [
    # GPS-like isotropic case, at LEO ≈ 20m and 0.02m/s
    [
        8.16e-12  0         0         0         0         0       ;
        0         8.16e-12  0         0         0         0       ;
        0         0         8.16e-12  0         0         0       ;
        0         0         0         8.16e-12  0         0       ;
        0         0         0         0         8.16e-12  0       ;
        0         0         0         0         0         8.16e-12;
    ],

    # Ground observation with moderate along-track uncertainty and T-v_r and R-v_t correlation
    [
        1e-10   0        0      0        -5e-11  0    ;
        0       2.5e-9   0      1.25e-9  0       0    ;
        0       0        1e-10  0        0       0    ;
        0       1.25e-9  0      2.5e-9   0       0    ;
        -5e-11  0        0      0        1e-10   0    ;
        0       0        0      0        0       1e-10;
    ],

    # Ground observation with higher along-track uncertainty and T-v_r and R-v_t correlation
    [
        1e-8    0        0      0        -9e-9   0   ;
        0       1e-6     0      9e-7     0       0   ;
        0       0        1e-8   0        0       0   ;
        0       9e-7     0      1e-6     0       0   ;
        -9e-9   0        0      0        1e-8    0   ;
        0       0        0      0        0       1e-8;
    ],
]
#! format: on

function case_covariance(μ, Σ)
    # Returns both Σ and the rotation matrix from RTN->ECI
    rv = mee_to_euclid(μ..., GM_EARTH)
    r = rv[1:3]
    v = rv[4:6]

    # RTN basis
    n = cross(r, v)
    t = cross(n, r)
    ru = r / norm(r)
    nu = n / norm(n)
    tu = t / norm(t)
    R_3x3 = hcat(ru, tu, nu)
    R = [R_3x3 zeros(3, 3); zeros(3, 3) R_3x3]

    # Scale matrix
    vorb = sqrt(GM_EARTH / norm(r))
    S = Diagonal([norm(r), norm(r), norm(r), vorb, vorb, vorb])

    # Bring the uncertainty to ECI
    return Symmetric(R * S * Σ * S' * R'), R
end

function main()
    GLMakie.activate!()
    for case_i = 1:length(CASES_MEAN_MEE)
        fig = Figure(size = (600, 600))
        ax3 = Axis3(
            fig[1, 1],
            title = "Case $case_i",
            xlabel = "x [km]",
            ylabel = "y [km]",
            zlabel = "z [km]",
            aspect = :data,
            xgridcolor = :gray,
            ygridcolor = :gray,
            zgridcolor = :gray,
        )


        μ = CASES_MEAN_MEE[case_i]
        Σ, R_eci2rtn = case_covariance(μ, CASES_COV_RTN_SCALED[case_i])

        @info "Case" case_i μ Σ

        plot_rv_ellipsoid(ax3, μ, Σ; relative = true, R = R_eci2rtn)

        # display(fig)
        display(GLMakie.Screen(), fig)
    end
end

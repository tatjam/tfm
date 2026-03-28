# MonteCarlo.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# Simple Monte-Carlo runner, given a starting point cloud, returns the
# samples of the system
# 

"""
   run_monte_carlo(p, samples, delta_t; reltol, abstol) 

Propagates the force model for each of the samples, returning the array of
ending positions after the given time has passed.
"""
function run_monte_carlo(
    p::ForceModel,
    samples::AbstractVector{<:SVector{6,T}},
    delta_t;
    reltol=1e-10,
    abstol=1e-10) where {T<:Real}

    n = length(samples)
    results = Vector{SVector{6,T}}(undef, n)

    Threads.@threads for i in eachindex(samples)
        results[i] = propagate_orbit(p, samples[i], delta_t, reltol=reltol, abstol=abstol)
    end

    return results
end

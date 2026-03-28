using Test
using StaticArrays
using DifferentialEquations

using OrbitalUncertainty

@testset "Dynamics" begin
    include("test_dynamics.jl")
end

@testset "Monte Carlo" begin
    include("test_monte_carlo.jl")
end


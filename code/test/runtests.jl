using Test
using StaticArrays
using DifferentialEquations

using OrbitalUncertainty

@testset "Dynamics" begin
    include("test_dynamics.jl")
end


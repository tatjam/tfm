using Test
using StaticArrays
using DifferentialEquations
using Distributions
using LinearAlgebra
using SatelliteToolbox

using OrbitalUncertainty

@testset "Dynamics" begin
    include("test_dynamics.jl")
end

@testset "Monte Carlo" begin
    include("test_monte_carlo.jl")
end

@testset "Unscented Transform" begin
    include("test_ut.jl")
end

@testset "GMV" begin
    include("test_gvm.jl")
end

@testset "GMV Distribution" begin
    include("test_gvm_dist.jl")
end

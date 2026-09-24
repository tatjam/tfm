using Random: MersenneTwister, shuffle
using LinearAlgebra: dot
using Statistics: cov


rng = MersenneTwister(1)
n, m, d = 300, 200, 6
X = randn(rng, d, n) .* [1.0, 5.0, 1.0, 0.1, 2.0, 1.0]
Y = randn(rng, d, m) .+ 0.3

@testset "Fast metric matches naive implementation" begin
    @test energy_metric(X, Y) ≈ energy_metric_reference(X, Y) rtol = 1e-8

    # a random relabelling, through the weight-vector identity
    D = distance_matrix(hcat(X, Y))
    p = shuffle(rng, 1:(n + m))
    idxI, idxJ = sort(p[1:n]), sort(p[(n + 1):end])
    ref = energy_metric_on_distance_matrix(D, idxI, idxJ, sum(D))
    w = fill(-1 / m, n + m)
    w[idxI] .= 1 / n
    @test -dot(w, center_distance_matrix!(copy(D)) * w) ≈ ref rtol = 1e-8
end

@testset "Float32 agrees with Float64" begin
    @test energy_metric(X, Y; precision = Float32) ≈ energy_metric(X, Y) rtol = 1e-3
end

@testset "Whitening" begin
    Xw, Yw = withen_union(X, Y)
    C = cov(hcat(Xw, Yw), dims = 2)
    @test maximum(abs, C .- (1:d .== (1:d)')) < 1e-8
    @test size(Xw) == size(X)          # non-mutating version
    @test X != Xw
end

@testset "Null distribution independent of batch size" begin
    a, b = zeros(300), zeros(300)
    energy_test_null_distribution!(a, X, Y; batch = 7, rng = MersenneTwister(2))
    energy_test_null_distribution!(b, X, Y; batch = 128, rng = MersenneTwister(2))
    @test a ≈ b rtol = 1e-8
end

@testset "P-values" begin
    Xa, Ya = randn(rng, d, n), randn(rng, d, m)
    @test energy_test(Xa, Ya, 500; rng = rng) > 0.01              # H0 true
    @test energy_test(Xa, Ya .+ 0.5, 500; rng = rng) == 1 / 501   # H0 false: min p
end


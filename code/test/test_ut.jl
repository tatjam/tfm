@testset "Sigma weights sum to one" begin
    α = 1e-3
    κ = 0.5
    β = 0.5


    @testset "Sane scaling parameters" begin
        α = 1e-3
        κ = 0
        β = 2
        sw = SigmaVectors(zeros(6), Diagonal(ones(6)), α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end

    @testset "Weird scaling parameters" begin
        sw = SigmaVectors(zeros(6), Diagonal(ones(6)), α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end

    @testset "Random mean" begin
        sw = SigmaVectors(rand(6), Diagonal(ones(6)), α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end

    @testset "Random covariance" begin
        # Note that covariance needs to be symmetric and positive definite 
        mat = rand(6, 6)
        P = mat * mat' + Diagonal(ones(6))
        sw = SigmaVectors(zeros(6), P, α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end

    @testset "From Distributions.jl" begin
        dist = MvNormal(zeros(6), I)
        sw = SigmaVectors(mean(dist), cov(dist), α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end

end
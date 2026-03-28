@testset "Sigma weights sum to one" begin
    @testset "Normal distribution, sane defaults" begin
        dist = MvNormal(zeros(6), I)
        α = 1e-3
        κ = 0
        β = 2
        sw = SigmaVectors(dist, α, κ, β)
        @test sum(sw.W) ≈ 1.0
    end
end
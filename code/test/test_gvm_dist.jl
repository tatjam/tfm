@testset "Basics" begin
    μ = [1.0, 2.0, 3.0]
    P = [1.0 0.0 0.0;
         0.0 1.0 0.0;
         0.0 0.0 1.0]
    α = 0.3
    β = [0.1, 0.2, 0.3]
    Γ = [0.0 0.0 0.0;
         0.0 0.0 0.0;
         0.0 0.0 0.0]
    κ = 0.1
    

    gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
    @testset "Constructor" begin 
        @test gvm.μ == μ
        @test gvm.α == α
        @test gvm.β == β
        @test gvm.Γ == Γ
        @test gvm.κ == κ

        @test istril(gvm.A)
        @test P ≈ gvm.A * gvm.A'
    end

    @testset "Sampling" begin
        vs = rand(gvm, 1000)
        @test all(l -> l == length(μ), length.(getindex.(vs, 1)))
        @test all(x -> all(isfinite, x), getindex.(vs, 1))
        @test all(θ -> -π <= θ <= π, getindex.(vs, 2))
    end
end

@testset "Marginals" begin
    μ = [1.0, 2.0, 3.0]
    P = [5.0 2.0 3.0;
         2.0 5.0 4.0;
         3.0 4.0 5.0]
    α = 0.3
    β = [1.5, 2.5, 3.5]
    Γ = [0.1 0.2 0.3;
         0.2 1.1 0.2;
         0.3 1.2 0.1]
    κ = 0.1

    n = 1000000

    @testset "Marginal on x is Gaussian" begin
        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
        xs = getindex.(rand(gvm, n), 1)
        @test mean(xs) ≈ μ atol=0.05
        @test cov(xs) ≈ P atol=0.05
    end

    @testset "Marginal on θ has angular average α if β and Γ are zero" begin
        gvm = GaussVonMises(μ, α, β * 0.0, Γ * 0.0, κ, P=P)
        θs = getindex.(rand(gvm, n), 2)
        avg = 1/n * sum(map(θ -> exp(im * θ), θs))
        avgangle = atan(imag(avg), real(avg)) 
        @test avgangle ≈ α atol=0.05
    end

    @testset "Marginal on θ is not Von Mises if β and Γ are non zero" begin
        # When β and Γ are present, the marginal is a "mixture" of Von Mises distributions,
        # which no longer is a Von Mises distribution
        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
        θs = getindex.(rand(gvm, n), 2)
        avg = 1/n * sum(map(θ -> exp(im * θ), θs))
        avgangle = atan(imag(avg), real(avg)) 
        @test abs(avgangle - α) > 0.1
    end
end


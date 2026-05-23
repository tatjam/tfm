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
        @test all(θ -> -π <= θ <= π, getindex.(vs, 4))
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
        xs = reduce(hcat, rand(gvm, n))[1:3, :]
        @test mean(xs,dims=2) ≈ μ atol=0.05
        @test cov(xs,dims=2) ≈ P atol=0.05
    end

    @testset "Marginal on θ has angular average α if β and Γ are zero" begin
        gvm = GaussVonMises(μ, α, β * 0.0, Γ * 0.0, κ, P=P)
        θs = reduce(hcat, rand(gvm, n))[4, :]
        avg = 1/n * sum(map(θ -> exp(im * θ), θs))
        avgangle = atan(imag(avg), real(avg)) 
        @test avgangle ≈ α atol=0.05
    end

    @testset "Marginal on θ is not Von Mises if β and Γ are non zero" begin
        # When β and Γ are present, the marginal is a "mixture" of Von Mises distributions,
        # which no longer is a Von Mises distribution
        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
        θs = reduce(hcat, rand(gvm, n))[4,:]
        avg = 1/n * sum(map(θ -> exp(im * θ), θs))
        avgangle = atan(imag(avg), real(avg)) 
        @test abs(avgangle - α) > 0.1
    end
end


@testset "Canonical Mahalanobis" begin
    μ = [0.0, 0.0, 0.0]
    P = [1.0 0.0 0.0;
         0.0 1.0 0.0;
         0.0 0.0 1.0]
    α = 0.0
    β = [0.0, 0.0, 0.0]
    Γ = [0.0 0.0 0.0;
         0.0 0.0 0.0;
         0.0 0.0 0.0]
    n = 100000

    @testset "Mahalanobis matches Canonical Mahalanobis for canonical GVM" begin
        κ = 1000.0
        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)

        x = [3.0; 4.0;5.0;1.0]

        @test mahalanobis(x, gvm) ≈ canon_mahalanobis(x[1:end-1], x[end], κ)
    end

    @testset "Expected Canonical Mahalanobis of samples with big κ" begin
        κ = 1000.0
        # If we average enough samples, their Mahalanobis distance will follow
        # 𝔼[zᵀz] + 𝔼[4κsin(½ϕ)] = length(μ) + 2κ(1 - I₁/I₀) ≈ 6 for big κ

        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
        samples = reduce(hcat, rand(gvm, n))
        avg_mahalanobis = mean(map(eachcol(samples)) do v
            mahalanobis(v, gvm)
        end)

        @test avg_mahalanobis ≈ 4 rtol=0.01
    end

    @testset "Canonical Mahalanobis of samples" begin
        κ = 0.5
        # Same as before, but we can't use the approximation for κ

        gvm = GaussVonMises(μ, α, β, Γ, κ, P=P)
        samples = reduce(hcat, rand(gvm, n))
        avg_mahalanobis = mean(map(eachcol(samples)) do v
            mahalanobis(v, gvm)
        end)

        ang_term = 2κ * b12(κ)[1]

        @test avg_mahalanobis ≈ 3 + ang_term rtol=0.01

    end

end

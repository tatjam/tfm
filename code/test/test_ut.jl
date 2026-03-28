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

@testset "Propagation" begin
    x0 = [1000.0e3, 2000.0e3, 3000.0e3, 1000.0, 5000.0, 10000.0]
    x0nospeed = [1000.0e3, 2000.0e3, 3000.0e3, 0.0, 0.0, 0.0]
    covmat = Diagonal(ones(6))

    # Velocity uncertainty can't be 0, otherwise the matrix doesn't have sqrt()
    covmatnospeed = Diagonal([ones(3); fill(1e-10, 3)])

    @testset "Empty force function alongside static point doesn't change the mean nor covariance" begin
        fm = ForceModel(())
        dist = run_ut(fm, x0nospeed, covmatnospeed, 1000)
        @test mean(dist) ≈ x0nospeed
        # Note, because we can't set velocity uncertainty to 0, we have some growth here
        @test cov(dist) ≈ covmatnospeed atol = 1e-3
    end

    @testset "Empty force function doesn't change the covariance nor speed" begin
        fm = ForceModel(())

        dist = run_ut(fm, x0, covmat, 1000)
        @test mean(dist)[4:end] ≈ x0[4:end]
        @test norm(cov(dist) - covmat) > 1.0
    end

    @testset "Keplerian propagation changes the distribution" begin
        fm = ForceModel(())

        dist = run_ut(fm, x0, covmat, 1000)
        @test norm(mean(dist) - x0) > 1000.0
        @test norm(cov(dist) - covmat) > 1.0
    end


    @testset "Very short orbital case matches sufficient Monte Carlo samples" begin
        Δt = 60
        nsamples = 100000

        # A sufficiently short Keplerian problem looks very linear
        dist = run_ut(EARTH_FM_WITH_J2, x0, covmat, Δt)

        # Monte-Carlo samples
        start_dist = MvNormal(x0, covmat)
        samples = rand(start_dist, nsamples)
        samples_svec = reinterpret(SVector{6,eltype(samples)}, vec(samples))
        samples_mc = run_monte_carlo(EARTH_FM_WITH_J2, samples_svec, Δt)

        # Flatten to a conventional matrix 
        samples_mc_mat = reduce(hcat, samples_mc)

        μ_mc = vec(mean(samples_mc_mat, dims=2))
        P_mc = cov(samples_mc_mat, dims=2)

        @test μ_mc ≈ mean(dist) rtol = 1e-6

        rel_error = norm(cov(dist) - P_mc) / norm(P_mc)
        @test rel_error < 0.01

    end
end
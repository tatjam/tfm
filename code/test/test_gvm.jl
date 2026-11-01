vm_ratio(κ, k) = besselix(k, κ) / besselix(0, κ)          # Iₖ(κ)/I₀(κ)

# Random lower-triangular factor with positive diagonal
function rand_chol(rng, n; scale = 0.3)
    M = scale .* randn(rng, n, n)
    return LowerTriangular(tril(M) + Diagonal(0.5 .+ rand(rng, n)))
end

function rand_gvm(rng, n; κ = 500.0, γscale = 1.0, with_Γ = true)
    A = rand_chol(rng, n)
    G = 0.3 .* randn(rng, n, n)
    Γ = with_Γ ? (G + G') / 2 : zeros(n, n)
    return GaussVonMises(randn(rng, n), randn(rng), γscale .* randn(rng, n), Γ, κ, A = A)
end

# rand(d, N) may return a matrix or a vector of vectors, handle both
cols(S::AbstractMatrix) = collect(eachcol(S))
cols(S::AbstractVector) = S

cond_mean(d, z) = d.α + dot(d.β, z) + dot(z, d.Γ, z) / 2

function angular_residual(d, x)
    z = d.A \ (x[1:(end-1)] .- d.μ)
    return x[end] - cond_mean(d, z)
end

# Covariance of the joint Gaussian (x, θ) that a GVM with Γ = 0 and large κ approximates
function joint_cov(A, β, σc2)
    Aβ = A * β
    return Symmetric([Matrix(A * A') Aβ; Aβ' (dot(β, β)+σc2)])
end

function matched(rng, n, σc2)
    A = rand_chol(rng, n)
    β = randn(rng, n)
    μe = randn(rng, n)
    α = 0.4
    Σ = joint_cov(A, β, σc2)
    g = moment_match(MvNormal([μe; α], Σ), n + 1)
    return (; g, A, β, μe, α, Σ)
end

function test_same_dist(out, ref; atol = 1e-8)
    @test out.μ ≈ ref.μ atol = atol
    @test out.α ≈ ref.α atol = atol
    @test Matrix(out.A) ≈ Matrix(ref.A) atol = atol
    @test out.β ≈ ref.β atol = atol
    @test Matrix(out.Γ) ≈ Matrix(ref.Γ) atol = atol
    @test out.κ ≈ ref.κ
end

unitvec(N, k) = SVector{N}(ntuple(i -> Float64(i == k), N))

@testset "canon_mahalanobis" begin
    z = [0.3, -1.2, 0.7]
    κ = 25.0
    @test canon_mahalanobis(zeros(3), 0.0, κ) == 0.0
    @test canon_mahalanobis(z, 0.0, κ) ≈ dot(z, z)
    @test canon_mahalanobis(z, 2π, κ) ≈ dot(z, z)               # periodic
    @test canon_mahalanobis(z, π, κ) ≈ dot(z, z) + 4κ           # antipode: chord = 2
    @test canon_mahalanobis(z, 0.4, κ) ≈ canon_mahalanobis(z, -0.4, κ)
    ϕ = 1e-3                                                     # small angle: κϕ²
    @test canon_mahalanobis(z, ϕ, κ) ≈ dot(z, z) + κ * ϕ^2 rtol = 1e-6
end

@testset "b12" begin
    for κ in (50.0, 1e3, 1e4)
        b1, b2 = b12(κ)
        @test b1 ≈ 1 - vm_ratio(κ, 1)
        @test b1 ≈ 1 / (2κ) rtol = 2 / κ                         # leading order
        # exact identity from I₀ - I₂ = (2/κ) I₁
        @test b2 ≈ (2 / κ) * (1 - b1) rtol = 1e-8
    end
end

@testset "mahalanobis" begin
    rng = MersenneTwister(1)
    d = rand_gvm(rng, 4)
    n = length(d.μ)

    @test mahalanobis([d.μ; d.α], d) ≈ 0 atol = 1e-12

    @testset "on the conditional-mean ridge d² = |z|² exactly" begin
        for _ = 1:20
            z = randn(rng, n)
            x = [d.μ .+ d.A * z; cond_mean(d, z)]
            @test mahalanobis(x, d) ≈ dot(z, z) rtol = 1e-10
        end
    end

    @testset "2π-periodic in the angle" begin
        x = [d.μ .+ 0.1 .* randn(rng, n); d.α + 0.2]
        x2 = copy(x)
        x2[end] += 2π
        x3 = copy(x)
        x3[end] -= 4π
        @test mahalanobis(x2, d) ≈ mahalanobis(x, d)
        @test mahalanobis(x3, d) ≈ mahalanobis(x, d)
    end

    @testset "non-negative" begin
        xs = [[d.μ .+ randn(rng, n); d.α + 3randn(rng)] for _ = 1:200]
        @test all(x -> mahalanobis(x, d) ≥ 0, xs)
    end

    @testset "large κ, Γ = 0: equals the joint Gaussian Mahalanobis" begin
        n = 3
        A = rand_chol(rng, n)
        β = randn(rng, n)
        κ = 1e6
        g = GaussVonMises(randn(rng, n), 0.7, β, zeros(n, n), κ, A = A)
        Σinv = inv(joint_cov(A, β, 1 / κ))
        m = [g.μ; g.α]
        for _ = 1:20
            z = randn(rng, n)
            x = [g.μ .+ A * z; g.α + dot(β, z) + 1e-3 * randn(rng)]
            δ = x .- m
            @test mahalanobis(x, g) ≈ dot(δ, Σinv * δ) rtol = 1e-5
        end
    end
end

@testset "sampling" begin
    Random.seed!(11)
    n, κ, N = 3, 100.0, 100_000
    d = rand_gvm(MersenneTwister(6), n; κ)
    S = cols(rand(d, N))

    # Euclidean marginal is N(μ, AAᵀ)
    X = reduce(hcat, [x[1:n] for x in S])
    @test vec(mean(X, dims = 2)) ≈ d.μ atol = 0.02
    @test cov(X, dims = 2) ≈ Matrix(d.A * d.A') atol = 0.03

    # Angular residual is von Mises(κ): E cos ϕ = I₁/I₀
    c = [cos(angular_residual(d, x)) for x in S]
    @test abs(mean(c) - vm_ratio(κ, 1)) < 6 * std(c) / sqrt(N)

    # Exact expectation of the distance under the model: n + 2κ(1 - I₁/I₀)
    d2 = [mahalanobis(x, d) for x in S]
    expected = n + 2κ * (1 - vm_ratio(κ, 1))
    @test abs(mean(d2) - expected) < 6 * std(d2) / sqrt(N)
end

@testset "GVMSigmaVectors κ = $κ" for κ in (5.0, 100.0, 1e4)
    rng = MersenneTwister(2)
    d = rand_gvm(rng, 3; κ)
    n = length(d.μ)
    σ = GVMSigmaVectors(d)
    W, X = σ.W, σ.χ
    @test size(X) == (n + 1, 2n + 3)

    idx = 1:size(X, 2)
    zs = [d.A \ (X[1:n, i] - d.μ) for i in idx]
    ϕs = [X[n+1, i] - cond_mean(d, zs[i]) for i in idx]

    @test sum(W) ≈ 1
    # Euclidean part reproduces N(0, I) in canonical coordinates, incl. 4th moment
    @test sum(W[i] * zs[i] for i in idx) ≈ zeros(n) atol = 1e-10
    @test sum(W[i] * zs[i] * zs[i]' for i in idx) ≈ I(n) atol = 1e-9
    @test [sum(W[i] * zs[i][j]^4 for i in idx) for j = 1:n] ≈ fill(3.0, n) atol = 1e-9
    # Angular part reproduces the first two circular moments of VM(κ)
    @test abs(sum(W[i] * sin(ϕs[i]) for i in idx)) < 1e-10
    @test sum(W[i] * cos(ϕs[i]) for i in idx) ≈ vm_ratio(κ, 1) rtol = 1e-8
    @test sum(W[i] * cos(2ϕs[i]) for i in idx) ≈ vm_ratio(κ, 2) rtol = 1e-8
end

# Very large κ is what you get for orbital states (angle nearly determined by the
# other elements). 1 - I₁/I₀ suffers catastrophic cancellation here, so these are
# expected to FAIL until the asymptotic expansion is used for large κ.
@testset "GVMSigmaVectors large κ = $κ" for κ in (1e5, 1e6, 1e8, 1e10, 1e12)
    d = rand_gvm(MersenneTwister(3), 3; κ)
    n = length(d.μ)
    σ = GVMSigmaVectors(d)
    @test all(isfinite, σ.χ)
    @test all(isfinite, σ.W)
    @test sum(σ.W) ≈ 1
    ϕs = [σ.χ[n+1, i] - cond_mean(d, d.A \ (σ.χ[1:n, i] - d.μ)) for i = 1:size(σ.χ, 2)]
    @test maximum(abs, ϕs) ≈ sqrt(3 / κ) rtol = 1e-2       # η → √(3/κ)
    @test σ.W[2] ≈ 1 / 6 rtol = 1e-2                       # w_η → 1/6
end

@testset "moment_match" begin
    rng = MersenneTwister(4)
    n = 4

    @testset "recovers structure, σc² = $σc2" for σc2 in (1.0, 1e-2, 1e-4, 1e-6)
        r = matched(rng, n, σc2)
        g = r.g
        @test g.μ ≈ r.μe
        @test g.α ≈ r.α
        @test Matrix(g.A) ≈ Matrix(r.A) rtol = 1e-8
        @test g.β ≈ r.β rtol = 1e-8
        @test all(iszero, g.Γ)
        # Defining equation I₁(κ)/I₀(κ) = exp(-σc²/2), tested on 1 - ratio
        # (the large-κ branch is an approximation, hence the loose rtol)
        @test (1 - vm_ratio(g.κ, 1)) ≈ -expm1(-σc2 / 2) rtol = 1e-3
        # Small conditional variance: κ ≈ 1/σc²
        σc2 ≤ 1e-2 && @test g.κ ≈ 1 / σc2 rtol = σc2
    end

    @testset "κ decreases as the conditional variance grows" begin
        κs = [matched(rng, n, s).g.κ for s in (1e-6, 1e-4, 1e-2, 1.0)]
        @test issorted(κs, rev = true)
    end

    @testset "no jump where the exact / asymptotic branches switch (0.999)" begin
        κs = map((0.999 - 1e-7, 0.999 + 1e-7)) do target
            matched(rng, n, -2 * log(target)).g.κ
        end
        @test κs[1] ≈ κs[2] rtol = 2e-3
    end

    @testset "GVM distance ≈ joint Gaussian distance for tiny σc²" begin
        σc2 = 1e-6
        r = matched(rng, n, σc2)
        Σinv = inv(r.Σ)
        m = [r.μe; r.α]
        for _ = 1:20
            z = randn(rng, n)
            x = [r.μe .+ r.A * z; r.α + dot(r.β, z) + sqrt(σc2) * randn(rng)]
            δ = x .- m
            @test mahalanobis(x, r.g) ≈ dot(δ, Σinv * δ) rtol = 1e-4
        end
    end

    @testset "degenerate conditional variance gives finite κ" begin
        r = matched(rng, n, 1e-10)
        # We have a limit at 1e8 due to special functions, but real value is higher
        @test isfinite(r.g.κ) && r.g.κ > 1e7
        @test isfinite(mahalanobis([r.μe; r.α], r.g))
        # sigma points must be usable at this κ (expected to fail until fixed)
        @test all(isfinite, GVMSigmaVectors(r.g).W)
    end
end

@testset "gvm_propagate with analytic maps" begin
    rng = MersenneTwister(5)
    n = 3
    d = rand_gvm(rng, n; κ = 300.0)

    @testset "identity" begin
        test_same_dist(gvm_propagate(identity, d), d)
    end

    @testset "angle shift" begin
        out = gvm_propagate(w -> w + 0.8 * unitvec(n + 1, n + 1), d)
        ref = GaussVonMises(d.μ, d.α + 0.8, d.β, d.Γ, d.κ, A = d.A)
        test_same_dist(out, ref)
    end

    @testset "angle scaling" begin
        # θ' = 2θ: a VM(κ) variable scaled by 2 is VM(κ/4), so the end
        # distribution is exact and κ must shrink by c² = 4, with
        # c = ∂f_θ/∂θ − β'ᵀA'⁻¹(∂f_x/∂θ) = 2 here.
        out = gvm_propagate(w -> SA[w[1], w[2], w[3], 2 * w[4]], d)
        ref = GaussVonMises(d.μ, 2 * d.α, 2 * d.β, 2 * d.Γ, d.κ / 4, A = d.A)
        test_same_dist(out, ref)
    end

    @testset "euclidean translation" begin
        c = randn(rng, n)
        sh = SVector{n + 1}(c..., 0.0)
        out = gvm_propagate(w -> w + sh, d)
        ref = GaussVonMises(d.μ + c, d.α, d.β, d.Γ, d.κ, A = d.A)
        test_same_dist(out, ref)
    end

    @testset "general linear map: exact closed-form result" begin
        # x' = Jx + c ;  θ' = θ + bᵀx + d0
        J = Matrix(1.0I, n, n) + 0.3 * randn(rng, n, n)
        c, b, d0 = randn(rng, n), randn(rng, n), 0.25
        Fm = SMatrix{n + 1,n + 1}([J zeros(n); b' 1.0])
        sh = SVector{n + 1}(c..., d0)
        out = gvm_propagate(w -> Fm * w + sh, d)

        A2 = cholesky(Symmetric(J * Matrix(d.A * d.A') * J')).L
        M = A2 \ (J * d.A)                       # maps z → z′, orthogonal here
        @test M * M' ≈ I(n) atol = 1e-10
        βo = d.β + d.A' * b
        ref = GaussVonMises(
            J * d.μ + c,
            d.α + d0 + dot(b, d.μ),
            M' \ βo,                             # M⁻ᵀ βo
            (M' \ d.Γ) / M,                      # M⁻ᵀ Γ M⁻¹
            d.κ,
            A = A2,
        )
        test_same_dist(out, ref; atol = 1e-7)
    end

    @testset "mildly nonlinear map is consistent with Monte Carlo" begin
        Random.seed!(21)
        n = 2
        A = LowerTriangular(Matrix(1.0I, n, n))
        d0 = GaussVonMises(zeros(n), 0.3, [0.5, -0.2], zeros(n, n), 200.0, A = A)
        f(w) = SA[w[1]+0.05*w[1]^2, w[2]+0.05*w[1]*w[2], w[3]+0.5*w[1]+0.02*w[2]^2]
        out = gvm_propagate(f, d0)
        Y = [f(w) for w in cols(rand(d0, 50_000))]
        d2 = [mahalanobis(y, out) for y in Y]
        @test mean(d2) ≈ n + 1 rtol = 0.05
        @test 0.90 < mean(d2 .<= quantile(Chisq(n + 1), 0.95)) < 0.99
    end
end

@testset "run_gvm over a very short step is ≈ identity" begin
    # Same imports / constants as your existing propagation test
    u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    mee = kepler_to_mee(euclid_to_kepler(u0..., GM_EARTH)...)
    d = GaussVonMises(mee[1:5], mee[6], zeros(5), zeros(5, 5), 100.0, P = I(5) * 0.001)
    # Pure two-body: the Euclidean MEE elements are exactly constant, so the
    # propagated distribution must equal the prior. J2 can't be used here: its
    # short-period rates have huge state gradients (e.g. ∂ṗ/∂k ≈ 56 m/s per
    # unit k at this state), so even Δt = 1e-3 s changes the covariance by
    # O(5-10%). The map is ≈ identity, the covariance is not.
    out = run_gvm(EARTH_FM_KEPLER, d, 1e-3)
    test_same_dist(out, d; atol = 1e-4)
end

@testset "Kepler AD derivatives match closed form" begin
    # Under pure Kepler the Euclidean MEE elements are constant and
    # L̇ = √(μp) w²/p² with w constant along the flow, so f_θ = L + Δt·L̇
    # up to O(Δt²), and (with p in place of the expert's a):
    #   ∂f_θ/∂p = −(3/2)(L̇/p)Δt,   ∂²f_θ/∂p² = (15/4)(L̇/p²)Δt.
    u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    mee = kepler_to_mee(euclid_to_kepler(u0..., GM_EARTH)...)
    Δt = 1e-3
    rate = force_model(mee, EARTH_FM_KEPLER, 0.0)[6]

    δf = ForwardDiff.jacobian(v -> propagate_orbit(EARTH_FM_KEPLER, v, Δt), mee)
    @test δf[6, 1] ≈ -(3 / 2) * (rate / mee[1]) * Δt rtol = 1e-6

    δ² = ForwardDiff.hessian(v -> propagate_orbit(EARTH_FM_KEPLER, v, Δt)[6], mee)
    @test δ²[1, 1] ≈ (15 / 4) * (rate / mee[1]^2) * Δt rtol = 1e-6
end

@testset "node slopes vs AD slopes (J2)" begin
    # Statistical linear regression on the propagated nodes vs AD at the mode.
    # The expert's criterion for dropping the Hessians: at a mild horizon the
    # two must agree to a few percent. μ and A are slope-independent, so they
    # match exactly; κ differs by the Stein-quadrature bias −1/(2κ) in the
    # sinφ·φ moment (~1% at κ = 100, vanishing for orbital κ).
    u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    mee = kepler_to_mee(euclid_to_kepler(u0..., GM_EARTH)...)
    d = GaussVonMises(mee[1:5], mee[6], zeros(5), zeros(5, 5), 100.0, P = I(5) * 0.001)

    out_ad = run_gvm(EARTH_FM_WITH_J2_KEPLER, d, 1e-3)
    out_nd = run_gvm(EARTH_FM_WITH_J2_KEPLER, d, 1e-3; slopes = :nodes)

    @test out_nd.μ == out_ad.μ
    @test Matrix(out_nd.A) == Matrix(out_ad.A)
    @test out_nd.β ≈ out_ad.β rtol = 0.1 atol = 1e-8
    @test Matrix(out_nd.Γ) ≈ Matrix(out_ad.Γ) rtol = 0.1 atol = 1e-8
    @test out_nd.κ ≈ out_ad.κ rtol = 0.1
end

# Algebra the propagator relies on: under z′ = M z the angular conditional mean
# α + β·z + ½zᵀΓz keeps its value iff β′ = M⁻ᵀβ and Γ′ = M⁻ᵀΓM⁻¹.
# Move this into a function in src/ and call it from gvm_propagate, then test that
# function: the version that uses Mβ and MΓMᵀ is only right when M is orthogonal.
@testset "reparametrisation of (β, Γ) under z′ = M z" begin
    rng = MersenneTwister(8)
    n = 4
    M = randn(rng, n, n) + 3I                  # deliberately not orthogonal
    β = randn(rng, n)
    G = randn(rng, n, n)
    Γ = (G + G') / 2
    β2 = M' \ β
    Γ2 = (M' \ Γ) / M
    for _ = 1:10
        z = randn(rng, n)
        z2 = M * z
        @test dot(β2, z2) + dot(z2, Γ2, z2) / 2 ≈ dot(β, z) + dot(z, Γ, z) / 2
    end
end

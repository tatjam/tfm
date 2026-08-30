
@testset "Normal squared is ChiSquared" begin
    # Under the Hermite polynomials, the normal is expressed
    # as just (0, 1, 0) and chi squared is (1,0,√2), noting that
    # the basis is (1, x, 1/√2 (x²-1)), for a normal input ξ
    # Galerkin projection should handle this exactly up to float precision!
    f(x) = x^2

    basis = hermite_basis(Float64, 3)
    proj = galerkin(f, basis)

    @test proj[1] ≈ 1 atol = 1e-9
    @test proj[2] ≈ 0 atol = 1e-9
    @test proj[3] ≈ sqrt(2) atol = 1e-9
end

@testset "Normal squared is ChiSquared (useless coefficients)" begin
    # Same as before but with a lot of expected 0 coefficients
    f(x) = x^2

    basis = hermite_basis(Float64, 20)
    proj = galerkin(f, basis)

    @test proj[1] ≈ 1 atol = 1e-9
    @test proj[2] ≈ 0 atol = 1e-9
    @test proj[3] ≈ sqrt(2) atol = 1e-9
    @test all(isapprox.(proj[4:end], 0, atol = 1e-9))

end

@testset "Basis self-projection" begin
    N = 20
    basis = hermite_basis(Float64, N)
    
    # Note that the quadrature is 
    @testset "Basis index $k" for k in 1:length(basis)
        p_k(x) = eval_basis(basis, x)[k]
        
        proj = galerkin(p_k, basis)
        
        # The projection of a unit vector of the basis on itself is 1 only
        # at its index, 0 otherwise
        expected = zeros(Float64, length(basis), 1)
        expected[k, 1] = 1.0

        @test proj ≈ expected atol = 1e-9
    end
end

@testset "Cubic polynomial expansion" begin
    basis = hermite_basis(Float64, 4)
    proj = galerkin(x -> x^3, basis)

    @test proj[1] ≈ 0.0 atol = 1e-9 
    @test proj[2] ≈ 3.0 atol = 1e-9 
    @test proj[3] ≈ 0.0 atol = 1e-9 
    @test proj[4] ≈ sqrt(6) atol = 1e-9
    @test proj[5] ≈ 0.0 atol = 1e-9
end


@testset "Multi-dimensional output projection" begin
    basis = hermite_basis(Float64, 3)
    f(x) = [2.5*x - 1.0, x^2]
    
    proj = galerkin(f, basis)
    
    @test size(proj) == (4, 2)
    @test proj[:, 1] ≈ [-1.0, 2.5, 0.0, 0.0] atol = 1e-9
    @test proj[:, 2] ≈ [1.0, 0.0, sqrt(2), 0.0] atol = 1e-9
end

@testset "Float32 Type Preservation" begin
    basis32 = hermite_basis(Float32, 3)
    proj32 = galerkin(x -> x^2, basis32)

    @test eltype(proj32) === Float32
    @test proj32[1] ≈ 1.0f0 atol = 1e-5
    @test proj32[3] ≈ sqrt(2.0f0) atol = 1e-5
end

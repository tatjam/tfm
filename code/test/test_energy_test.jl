@testset "Test of identical univariate samples is 1" begin
    X = randn(10000)
    X = reshape(X, (1, length(X)))
    Y = copy(X)

    @test energy_test(X, Y, 1000) ≈ 1

end

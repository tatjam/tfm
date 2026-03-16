function propagate_orbit(fm::ForceModel, u0, t)
    tspan = (0.0, t)
    prob = ODEProblem(newton_model, u0, tspan, fm)
    sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10)
    SA[sol[1,end], sol[2,end], sol[3,end], sol[4,end], sol[5,end], sol[6,end]]
end

@testset "J2 compare against poliastro" begin
    μ = 398600441800000.0
    j2 = 0.00108263
    r = 6378136.6

    cbf = TwoBodyForce(μ)
    j2f = J2Force(μ, r, j2)
    fm = ForceModel((cbf, j2f))

    orbit1_u0 = SA[6771.358863, 0, 0, 0, 4.76807358, 6.01581168] .* 1e3
    orbit1_u1 = SA[6679.99961672, -711.28165725, -849.9863342, 1.25649055, 4.69949408, 5.93788535] .* 1e3

    @test propagate_orbit(fm, orbit1_u0, 90.0 * 60.0) ≈ orbit1_u1 rtol=1e-3

end

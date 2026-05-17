using Test
using LinearAlgebra
using SwitchingControl

# ============================================================
# SIMPLE TEST PROBLEM
# ============================================================

function f_before!(dy, y, p, t)
    x, A = y
    dy[1] = 1.0
    dy[2] = -0.5 * A
    return nothing
end

function f_after!(dy, y, p, t)
    x, A = y
    dy[1] = -0.25
    dy[2] = 0.1 * A
    return nothing
end

function switching_condition(y, p, t)
    return y[1] - 1.0
end

bc = BoundaryConditions([0.0, 10.0], nothing)

prob = SwitchingProblem(
    f_before!,
    f_after!,
    switching_condition,
    bc,
    nothing,
    (0.0, 5.0),
    1
)

@testset "SwitchingControl IVP Tests" begin

    sol = SwitchingControl.solve(prob; mode = :ivp)

    @test sol.sol.converged == true
    @test sol.sol.switched == true
    @test !isnothing(sol.sol.tau)

    @test sol.sol.solution_before !== nothing
    @test sol.sol.solution_after !== nothing

    dims_pre = unique(length(u) for u in sol.sol.solution_before.u)
    dims_post = unique(length(u) for u in sol.sol.solution_after.u)

    @test dims_pre == [2]
    @test dims_post == [2]

    u_before = sol.sol.solution_before.u[end]
    u_after  = sol.sol.solution_after.u[1]

    @test norm(u_after .- u_before) < 1e-12
    @test abs(u_before[1] - 1.0) < 1e-6

    @test all(isfinite, sol.sol.solution_after.u[end])

    @test length(sol.sol.solution_before.u) > 1
    @test length(sol.sol.solution_after.u) > 1
end
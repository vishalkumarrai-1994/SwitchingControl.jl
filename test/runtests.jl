using Test
using LinearAlgebra
using SwitchingControl

# ============================================================
# SIMPLE TEST PROBLEM
# ============================================================

#
# State vector:
#
# y[1] = x
# y[2] = A
#
# Switch occurs when x = 1
#

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

    x = y[1]

    return x - 1.0
end

# ============================================================
# BUILD PROBLEM
# ============================================================

bc = BoundaryConditions(

    [0.0, 10.0],   # initial condition
    nothing
)

prob = SwitchingProblem(

    f_before!,
    f_after!,
    switching_condition,
    bc,
    nothing,
    (0.0, 5.0),
    1
)

# ============================================================
# TEST SET
# ============================================================

@testset "SwitchingControl IVP Tests" begin

    sol = SwitchingControl.solve(
        prob;
        mode = :ivp
    )

    # --------------------------------------------------------
    # BASIC FLAGS
    # --------------------------------------------------------

    @test sol.sol.converged == true
    @test sol.sol.switched == true

    # --------------------------------------------------------
    # SWITCH TIME EXISTS
    # --------------------------------------------------------

    @test !isnothing(sol.sol.tau)

    # --------------------------------------------------------
    # PRE / POST SOLUTIONS EXIST
    # --------------------------------------------------------

    @test sol.sol.solution_before !== nothing
    @test sol.sol.solution_after !== nothing

    # --------------------------------------------------------
    # DIMENSION CONSISTENCY
    # --------------------------------------------------------

    dims_pre = unique(
        length(u) for u in sol.sol.solution_before.u
    )

    dims_post = unique(
        length(u) for u in sol.sol.solution_after.u
    )

    @test dims_pre == [2]
    @test dims_post == [2]

    # --------------------------------------------------------
    # CONTINUITY AT SWITCH
    # --------------------------------------------------------

    u_before = sol.sol.solution_before.u[end]
    u_after  = sol.sol.solution_after.u[1]

    jump_norm = norm(u_after .- u_before)

    @test jump_norm < 1e-12

    # --------------------------------------------------------
    # SWITCH CONDITION SATISFIED
    # --------------------------------------------------------

    x_switch = u_before[1]

    @test abs(x_switch - 1.0) < 1e-6

    # --------------------------------------------------------
    # TERMINAL STATE FINITE
    # --------------------------------------------------------

    uT = sol.sol.solution_after.u[end]

    @test all(isfinite, uT)

    # --------------------------------------------------------
    # TRAJECTORY NONEMPTY
    # --------------------------------------------------------

    @test length(sol.sol.solution_before.u) > 1
    @test length(sol.sol.solution_after.u) > 1

end

println("\n====================================")
println("All SwitchingControl tests passed.")
println("====================================")
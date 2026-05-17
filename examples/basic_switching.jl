using SwitchingControl
using CairoMakie

# ============================================================
# SIMPLE SWITCHING CONTROL EXAMPLE
#
# A particle moves with:
#
#   ẋ = v
#   v̇ = a
#
# Before switching:
#   a = +1
#
# After switching:
#   a = -2
#
# Switching occurs when:
#   x(t) = 10
#
# ============================================================

# ============================================================
# PRE-SWITCH DYNAMICS
# ============================================================

function dynamics_before!(dy, y, p, t)

    x, v = y

    dy[1] = v
    dy[2] = 1.0

    return nothing
end

# ============================================================
# POST-SWITCH DYNAMICS
# ============================================================

function dynamics_after!(dy, y, p, t)

    x, v = y

    dy[1] = v
    dy[2] = -2.0

    return nothing
end

# ============================================================
# SWITCHING CONDITION
# ============================================================

function switching_condition(y, p, t)

    x, v = y

    # switch when x = 10
    return x - 10.0
end

# ============================================================
# BOUNDARY CONDITIONS
# ============================================================

u0 = [
    0.0,   # initial position
    0.0    # initial velocity
]

bc = BoundaryConditions(u0, nothing)

# ============================================================
# BUILD PROBLEM
# ============================================================

prob = SwitchingProblem(
    dynamics_before!,
    dynamics_after!,
    switching_condition,
    bc,
    nothing,
    (0.0, 10.0),
    1
)

# ============================================================
# SOLVE
# ============================================================

sol = SwitchingControl.solve(prob; mode = :ivp)

# ============================================================
# BASIC OUTPUT
# ============================================================

println("\n====================================")
println("SwitchingControl Example")
println("====================================")

println("Switched?  : ", sol.sol.switched)
println("Switch time: ", sol.sol.tau)

# ============================================================
# EXTRACT TRAJECTORIES
# ============================================================

t_pre = sol.sol.solution_before.t
x_pre = [u[1] for u in sol.sol.solution_before.u]
v_pre = [u[2] for u in sol.sol.solution_before.u]

if sol.sol.switched

    t_post = sol.sol.solution_after.t
    x_post = [u[1] for u in sol.sol.solution_after.u]
    v_post = [u[2] for u in sol.sol.solution_after.u]

else

    t_post = Float64[]
    x_post = Float64[]
    v_post = Float64[]
end

# ============================================================
# PLOTS
# ============================================================

fig = Figure(size = (900, 400))

# ------------------------------------------------------------
# POSITION
# ------------------------------------------------------------

ax1 = Axis(
    fig[1, 1],
    title = "Position",
    xlabel = "Time",
    ylabel = "x(t)"
)

lines!(ax1, t_pre, x_pre, linewidth = 3)

if sol.sol.switched
    lines!(ax1, t_post, x_post, linewidth = 3)
    vlines!(ax1, [sol.sol.tau], linestyle = :dash)
end

# ------------------------------------------------------------
# VELOCITY
# ------------------------------------------------------------

ax2 = Axis(
    fig[1, 2],
    title = "Velocity",
    xlabel = "Time",
    ylabel = "v(t)"
)

lines!(ax2, t_pre, v_pre, linewidth = 3)

if sol.sol.switched
    lines!(ax2, t_post, v_post, linewidth = 3)
    vlines!(ax2, [sol.sol.tau], linestyle = :dash)
end

display(fig)
using DifferentialEquations
using SciMLBase
using NonlinearSolve
using FiniteDiff

using .Types
using .Solution

# ============================================================
# PUBLIC DISPATCHER
# ============================================================
"""
Solve a switching dynamical system.

Supports IVP-based endogenous switching and experimental shooting methods.
"""
function solve(prob::SwitchingProblem; mode::Symbol = :ivp, kwargs...)

    if mode == :ivp

        return solve_ivp(prob)

    elseif mode == :shooting

        @assert haskey(kwargs, :guess) "Shooting mode requires guess = ..."

        return solve_shooting(prob, kwargs[:guess])

    elseif mode == :free_boundary

        @assert haskey(kwargs, :guess)     "free_boundary mode requires guess = ..."
        @assert haskey(kwargs, :make_u0)   "free_boundary mode requires make_u0 = ..."
        @assert haskey(kwargs, :residual!) "free_boundary mode requires residual! = ..."

        return solve_free_boundary(
            prob, kwargs[:guess];
            make_u0 = kwargs[:make_u0],
            residual! = kwargs[:residual!],
            n_residuals = get(kwargs, :n_residuals, length(kwargs[:guess]))
        )

    else

        error("Unknown mode: $mode. Use :ivp, :shooting, or :free_boundary.")

    end
end


# ============================================================
# INTERNAL HELPERS
# ============================================================

"""
    terminal_state(sol)

Extract terminal model state from a SwitchingSolutionView.
"""
function terminal_state(sol::SwitchingSolutionView)

    traj = sol.sol

    if traj.switched

        @assert traj.solution_after !== nothing

        return traj.solution_after.u[end]

    else

        return traj.solution_before.u[end]

    end
end


# ============================================================
# IVP + SWITCHING SOLVER
# ============================================================

function solve_ivp(prob::SwitchingProblem)

    @assert prob.max_switches == 1

    f1! = prob.dynamics_before!
    f2! = prob.dynamics_after!

    g = prob.switching_condition

    p = prob.parameters

    t0, Tf = prob.tspan

    bc = prob.boundary_conditions

    # --------------------------------------------------------
    # INITIAL STATE
    # --------------------------------------------------------

    u0 = copy(bc.initial)

    T = eltype(u0)

    switched = Ref(false)

    τ = Ref(T(NaN))

    u_switch = Ref(copy(u0))

    # --------------------------------------------------------
    # CALLBACK
    # --------------------------------------------------------

    function condition(u, t, integrator)

        return g(u, p, t)

    end

    function affect!(integrator)

        switched[] = true

        τ[] = integrator.t

        # preserve AD types
        u_switch[] = copy(integrator.u)

        terminate!(integrator)

        return nothing
    end

    cb = ContinuousCallback(condition, affect!)

    # --------------------------------------------------------
    # PRE-SWITCH SYSTEM
    # --------------------------------------------------------

    prob1 = ODEProblem(
        f1!,
        u0,
        (t0, Tf),
        p
    )

    sol1 = DifferentialEquations.solve(
        prob1,
        Tsit5(),
        callback = cb
    )

    # --------------------------------------------------------
    # NO SWITCH CASE
    # --------------------------------------------------------

    if !switched[]

        raw = SwitchingSolution(
            nothing,
            sol1,
            nothing,
            false,
            true,
            0.0
        )

        return SwitchingSolutionView(raw)

    end

    # --------------------------------------------------------
    # POST-SWITCH SYSTEM
    # --------------------------------------------------------

    prob2 = ODEProblem(
        f2!,
        u_switch[],
        (τ[], Tf),
        p
    )

    sol2 = DifferentialEquations.solve(
        prob2,
        Tsit5()
    )

    raw = SwitchingSolution(
        τ[],
        sol1,
        sol2,
        true,
        true,
        0.0
    )

    return SwitchingSolutionView(raw)
end


# ============================================================
# GENERAL FREE-BOUNDARY / SHOOTING SOLVER
# ============================================================

"""
    solve_free_boundary(prob, guess; make_u0, residual!, n_residuals = length(guess))

Solve for free parameters `θ` such that a user-supplied residual
vanishes, where the residual is evaluated on the switching solution
obtained by solving `prob`'s IVP dynamics with initial condition
`make_u0(θ)`.

This is deliberately more general than a classic two-point-BVP
shooting method:

- `make_u0(θ) -> u0` builds the initial state from the free
  parameters. `θ` need not be raw state components — it can be
  coefficients of a local series expansion near a singular point,
  a subset of state entries, or anything else that determines `u0`.

- `residual!(F, sol, θ, p)` fills the residual vector `F` from the
  full `SwitchingSolutionView` `sol`. This can reference the
  terminal state (recovering classic BVP shooting), or the state
  *at the switch itself* via `sol.sol.tau` and
  `sol.sol.solution_before.u[end]` — which is what endogenous
  free-boundary conditions such as value-matching / smooth-pasting
  require. The old behaviour (match a fixed terminal state at a
  fixed final time) is one instance of this; see `solve_shooting`
  below.

Returns the `NonlinearSolve` solution object; the converged
parameters are in `.u`, and calling `make_u0` on `.u` reconstructs
the boundary condition that solves the problem.
"""
function solve_free_boundary(prob::SwitchingProblem, guess;
                              make_u0,
                              residual!,
                              n_residuals::Int = length(guess))

    guess = collect(guess)

    function F!(F, θ, p)

        @assert θ !== nothing

        T = eltype(θ)

        u0 = T.(make_u0(θ))

        new_bc = BoundaryConditions(
            u0,
            prob.boundary_conditions.terminal
        )

        new_prob = SwitchingProblem(
            prob.dynamics_before!,
            prob.dynamics_after!,
            prob.switching_condition,
            new_bc,
            prob.parameters,
            prob.tspan,
            prob.max_switches
        )

        sol = solve_ivp(new_prob)

        residual!(F, sol, θ, prob.parameters)

        return nothing
    end

    F0 = zeros(eltype(guess), n_residuals)

    f = NonlinearFunction(
        F!;
        resid_prototype = F0
    )

    prob_nl = NonlinearProblem(
        f,
        guess,
        nothing
    )

    return NonlinearSolve.solve(
        prob_nl,
        NewtonRaphson(; autodiff = AutoFiniteDiff())
    )
end


# ============================================================
# CLASSIC TWO-POINT BVP SHOOTING (special case of the above)
# ============================================================

"""
Shoots on the trailing components of the initial state to match a
fixed terminal state at the problem's fixed final time. Kept for
backward compatibility; new problems — especially endogenous
free-boundary problems like smooth-pasting/value-matching — should
generally use `solve_free_boundary` directly, since it also lets
the residual be evaluated at the switching point rather than only
at the terminal time.
"""
function solve_shooting(prob::SwitchingProblem, guess)

    bc = prob.boundary_conditions

    @assert bc.terminal !== nothing \
        "Terminal condition required for shooting mode"

    xT = collect(bc.terminal)
    x0_fixed = collect(bc.initial)
    n = length(xT)

    make_u0 = θ -> begin
        T = eltype(θ)
        x0 = T.(x0_fixed)
        nf = length(θ)
        x0[end - nf + 1:end] .= θ
        x0
    end

    residual! = (F, sol, θ, p) -> begin
        T = eltype(θ)
        x_end_full = terminal_state(sol)
        x_end = x_end_full[1:n]
        F .= x_end .- T.(xT)
        return nothing
    end

    return solve_free_boundary(
        prob, guess;
        make_u0 = make_u0,
        residual! = residual!,
        n_residuals = n
    )
end
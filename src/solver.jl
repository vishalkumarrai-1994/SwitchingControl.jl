using DifferentialEquations
using SciMLBase
using NonlinearSolve
using FiniteDiff

using .Types
using .Solution

# ============================================================
# PUBLIC DISPATCHER
# ============================================================

function solve(prob::SwitchingProblem; mode::Symbol = :ivp, kwargs...)

    if mode == :ivp

        return solve_ivp(prob)

    elseif mode == :shooting

        @assert haskey(kwargs, :guess) "Shooting mode requires guess = ..."

        return solve_shooting(prob, kwargs[:guess])

    else

        error("Unknown mode: $mode. Use :ivp or :shooting.")

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
# SHOOTING METHOD
# ============================================================

function solve_shooting(prob::SwitchingProblem, guess)

    bc = prob.boundary_conditions

    @assert bc.terminal !== nothing \
        "Terminal condition required for shooting mode"

    xT = collect(bc.terminal)

    x0_fixed = collect(bc.initial)

    guess = collect(guess)

    n = length(xT)

    # =========================================================
    # RESIDUAL FUNCTION
    # =========================================================

    function residual!(F, x_free, p)

        @assert x_free !== nothing

        T = eltype(x_free)

        # ----------------------------------------------------
        # RECONSTRUCT INITIAL STATE
        # ----------------------------------------------------

        x0 = T.(x0_fixed)

        # assumes final entries are free variables
        nf = length(x_free)

        x0[end - nf + 1:end] .= x_free

        xT_local = T.(xT)

        new_bc = BoundaryConditions(
            x0,
            xT_local
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

        x_end_full = terminal_state(sol)

        # compare only model dimensions
        x_end = x_end_full[1:n]

        F .= x_end .- xT_local

        return nothing
    end

    # =========================================================
    # NONLINEAR PROBLEM
    # =========================================================

    F0 = zeros(eltype(guess), n)

    f = NonlinearFunction(
        residual!;
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
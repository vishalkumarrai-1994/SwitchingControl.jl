using DifferentialEquations
using SciMLBase
using NonlinearSolve
using FiniteDiff

using .Types
using .Solution

# ============================================================
# ALGORITHM RESOLUTION (lets callers pass plain Symbols instead of
# needing `using NonlinearSolve` themselves to name AutoFiniteDiff,
# TrustRegion, etc. Advanced users can still pass real NonlinearSolve
# / ADTypes objects directly if they've imported those packages --
# anything that isn't a recognized Symbol is passed through as-is.)
# ============================================================

function _resolve_autodiff(autodiff)

    autodiff isa Symbol || return autodiff

    if autodiff === :central
        return AutoFiniteDiff(fdtype = Val(:central))
    elseif autodiff === :forward_diff
        return AutoFiniteDiff(fdtype = Val(:forward))
    elseif autodiff === :forward
        return AutoForwardDiff()
    else
        error("Unknown autodiff = :$autodiff. Use :central, :forward_diff, " *
              "or :forward, or pass an ADTypes object directly.")
    end
end

function _resolve_nl_alg(nl_alg, autodiff)

    nl_alg isa Symbol || return nl_alg

    if nl_alg === :trust_region
        return TrustRegion(; autodiff = autodiff)
    elseif nl_alg === :newton
        return NewtonRaphson(; autodiff = autodiff)
    elseif nl_alg === :levenberg_marquardt
        return LevenbergMarquardt(; autodiff = autodiff)
    else
        error("Unknown nl_alg = :$nl_alg. Use :trust_region, :newton, " *
              ":levenberg_marquardt, or pass a NonlinearSolve algorithm directly.")
    end
end

# ============================================================
# LEG SPEC (must be defined before anything below uses it as a
# type annotation, e.g. solve(::Vector{Leg}, ...) and
# solve_multi_leg(::Vector{Leg}, ...))
# ============================================================

"""
    Leg(; dynamics!, tspan, parameters, switching_condition, make_u0, max_switches = 1)

Specification for one independently-integrated leg of a multi-leg
free-boundary problem — e.g. Dixit-style entry/exit, where separate
ODEs are integrated forward and backward over disjoint/overlapping
ranges of the state variable and glued together via matching
conditions at endogenous thresholds, rather than a single IVP with
one switch.

Fields:
- `dynamics!(du, u, p, t)`: right-hand side for this leg. The same
  function is used before *and* after any switch, since for these
  problems the `ContinuousCallback` exists purely to stop the
  integration at a threshold — it doesn't change the dynamics.
- `tspan::Tuple`: integration range for this leg. Can be reversed
  (e.g. `(P_far, P_lower)`) to integrate backward.
- `parameters`: model parameters `p`, passed through to `dynamics!`.
- `switching_condition(θ) -> (u, p, t) -> value`: given the current
  free-parameter vector `θ`, returns the scalar switching-condition
  function whose zero crossing stops this leg (e.g.
  `θ -> (u, p, t) -> t - θ[3]` to stop at `θ[3]`).
- `make_u0(θ) -> u0`: given `θ`, returns this leg's initial state.
- `max_switches::Int = 1`.
"""
struct Leg
    dynamics!::Function
    tspan::Tuple
    parameters::Any
    switching_condition::Function
    make_u0::Function
    max_switches::Int
end

Leg(; dynamics!, tspan, parameters, switching_condition, make_u0, max_switches::Int = 1) =
    Leg(dynamics!, tspan, parameters, switching_condition, make_u0, max_switches)

# ============================================================
# PUBLIC DISPATCHER
# ============================================================
"""
Solve a switching dynamical system.

Supports IVP-based endogenous switching and experimental shooting methods.
"""
function solve(prob::SwitchingProblem; mode::Symbol = :ivp, kwargs...)

    if mode == :ivp

        return solve_ivp(prob; kwargs...)

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
            n_residuals = get(kwargs, :n_residuals, length(kwargs[:guess])),
            autodiff = get(kwargs, :autodiff, AutoFiniteDiff())
        )

    elseif mode == :multi_leg

        @assert haskey(kwargs, :legs)      "multi_leg mode requires legs = [...]"
        @assert haskey(kwargs, :guess)     "multi_leg mode requires guess = ..."
        @assert haskey(kwargs, :residual!) "multi_leg mode requires residual! = ..."

        return solve_multi_leg(
            kwargs[:legs], kwargs[:guess];
            residual! = kwargs[:residual!],
            n_residuals = get(kwargs, :n_residuals, length(kwargs[:guess])),
            autodiff = get(kwargs, :autodiff, :central),
            nl_alg = get(kwargs, :nl_alg, :trust_region),
            ivp_kwargs = get(kwargs, :ivp_kwargs, (;)),
            penalty = get(kwargs, :penalty, 1e6),
            nlsolve_kwargs = get(kwargs, :nlsolve_kwargs, (;))
        )

    else

        error("Unknown mode: $mode. Use :ivp, :shooting, :free_boundary, or :multi_leg.")

    end
end

"""
    solve(legs::Vector{Leg}, guess; residual!, kwargs...)

Convenience entry point for multi-leg free-boundary problems (see
`solve_multi_leg`) that doesn't require constructing a placeholder
`SwitchingProblem` just to call the main dispatcher. Equivalent to
`solve(dummy_prob; mode = :multi_leg, legs = legs, guess = guess,
residual! = residual!, kwargs...)`.
"""
function solve(legs::Vector{Leg}, guess;
                residual!,
                n_residuals::Int = length(guess),
                autodiff = :central,
                nl_alg = :trust_region,
                ivp_kwargs = (;),
                penalty = 1e6,
                nlsolve_kwargs = (;))

    return solve_multi_leg(
        legs, guess;
        residual! = residual!,
        n_residuals = n_residuals,
        autodiff = autodiff,
        nl_alg = nl_alg,
        ivp_kwargs = ivp_kwargs,
        penalty = penalty,
        nlsolve_kwargs = nlsolve_kwargs
    )
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

function solve_ivp(prob::SwitchingProblem; kwargs...)

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
        callback = cb;
        kwargs...
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
        Tsit5();
        kwargs...
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
    solve_free_boundary(prob, guess; make_u0, residual!, n_residuals = length(guess),
                         autodiff = AutoFiniteDiff(), nl_alg = NewtonRaphson(; autodiff),
                         ivp_kwargs = (;))

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

`autodiff` controls the differencing scheme used for the outer
NonlinearSolve Jacobian (default `AutoFiniteDiff()`); pass e.g.
`AutoFiniteDiff(fdtype = Val(:central))` or `AutoForwardDiff()` if
`make_u0`/`residual!`/the dynamics are dual-number safe.
`nl_alg` overrides the nonlinear solver entirely (e.g. `TrustRegion`
or `LevenbergMarquardt` when Newton stalls or diverges).
`ivp_kwargs` are forwarded to the inner `solve_ivp` call (e.g.
`(abstol = 1e-12, reltol = 1e-12)`) to keep ODE integration error
well below the outer Jacobian's finite-difference step size.

Returns the `NonlinearSolve` solution object; the converged
parameters are in `.u`, and calling `make_u0` on `.u` reconstructs
the boundary condition that solves the problem.
"""
function solve_free_boundary(prob::SwitchingProblem, guess;
                              make_u0,
                              residual!,
                              n_residuals::Int = length(guess),
                              autodiff = AutoFiniteDiff(),
                              nl_alg = NewtonRaphson(; autodiff = autodiff),
                              ivp_kwargs = (;))

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

        sol = solve_ivp(new_prob; ivp_kwargs...)

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
        nl_alg
    )
end


# ============================================================
# MULTI-LEG FREE-BOUNDARY SOLVER
# ============================================================

"""
    solve_multi_leg(legs::Vector{Leg}, guess;
                     residual!, n_residuals = length(guess),
                     autodiff = AutoFiniteDiff(), nl_alg = nothing,
                     ivp_kwargs = (;), penalty = 1e6)

Solve a joint free-boundary problem defined over several
independently-integrated `Leg`s, sharing one free-parameter vector
`θ`. Each call of the residual re-solves *every* leg's IVP (via
`solve_ivp`, so `dynamics!` need not be dual-number-safe unless you
pass `AutoForwardDiff()`), then hands the resulting
`Vector{SwitchingSolutionView}` — in the same order as `legs` — to:

    residual!(F, sols, θ, p)

`p` here is `legs[1].parameters` (legs are expected to share
parameters; pass whatever composite object you need through that
single field if they don't literally share one).

If any leg fails to switch for a given `θ` (i.e. its target
threshold lies outside the leg's integrated range), `F` is filled
with `penalty` and `residual!` is *not* called for that evaluation
— this keeps bad/out-of-range Newton or line-search iterates from
throwing, letting globalized algorithms (`TrustRegion`,
`LevenbergMarquardt`) backtrack instead of erroring out.

`nl_alg` defaults to `:trust_region` — globalization is the sensible
default for this problem class, since a plain Newton full step is
what tends to blow up into the `1e14`-sized iterates this kind of
matching system produces when given a rough guess. Other built-in
shorthands: `:newton`, `:levenberg_marquardt`. `autodiff` defaults
to `:central` (central finite differences); other shorthands:
`:forward_diff`, `:forward` (forward-mode AD, only valid if
`dynamics!`/`make_u0`/`residual!` are dual-number safe). None of
this requires `using NonlinearSolve` in caller code — the Symbols
are resolved internally. Advanced users who *have* imported
NonlinearSolve/ADTypes themselves can still pass real algorithm/AD
objects directly instead of a Symbol (e.g. `nl_alg =
LevenbergMarquardt(; autodiff = AutoForwardDiff())`); anything that
isn't a recognized Symbol is passed through unchanged.

`ivp_kwargs` (e.g. `(abstol = 1e-12, reltol = 1e-12)`) are forwarded
to every `solve_ivp` call — keep these tight relative to
`autodiff`'s finite-difference step, or the outer Jacobian will be
dominated by ODE integration noise rather than true sensitivity.
`nlsolve_kwargs` (e.g. `(abstol = 1e-9, maxiters = 100)`) are
forwarded to the outer `NonlinearSolve.solve` call itself.

Returns the `NonlinearSolve` solution object; `.u` holds the
converged `θ`.
"""
function solve_multi_leg(legs::Vector{Leg}, guess;
                          residual!,
                          n_residuals::Int = length(guess),
                          autodiff = :central,
                          nl_alg = :trust_region,
                          ivp_kwargs = (;),
                          penalty = 1e6,
                          nlsolve_kwargs = (;))

    guess = collect(guess)

    autodiff_resolved = _resolve_autodiff(autodiff)
    nl_alg_resolved    = _resolve_nl_alg(nl_alg, autodiff_resolved)

    function F!(F, θ, _p)

        T = eltype(θ)

        sols = Vector{SwitchingSolutionView}(undef, length(legs))

        for (i, leg) in enumerate(legs)

            u0 = T.(leg.make_u0(θ))
            bc = BoundaryConditions(u0, nothing)
            g  = leg.switching_condition(θ)

            leg_prob = SwitchingProblem(
                leg.dynamics!,
                leg.dynamics!,
                g,
                bc,
                leg.parameters,
                leg.tspan,
                leg.max_switches
            )

            sol = solve_ivp(leg_prob; ivp_kwargs...)

            if !sol.sol.switched
                fill!(F, T(penalty))
                return nothing
            end

            sols[i] = sol
        end

        residual!(F, sols, θ, legs[1].parameters)

        return nothing
    end

    F0 = zeros(eltype(guess), n_residuals)
    f = NonlinearFunction(F!; resid_prototype = F0)
    prob_nl = NonlinearProblem(f, guess, nothing)

    return NonlinearSolve.solve(prob_nl, nl_alg_resolved; nlsolve_kwargs...)
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
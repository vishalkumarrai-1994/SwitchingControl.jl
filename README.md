# SwitchingControl.jl

![CI](https://github.com/vishalkumarrai-1994/SwitchingControl.jl/actions/workflows/ci.yml/badge.svg)


A Julia framework for solving dynamic optimization and switching-control problems through event-driven numerical integration.

`SwitchingControl.jl` is designed for optimal control, economics, engineering, and hybrid dynamical systems where trajectories evolve under one regime until a switching condition is triggered, after which the system transitions to a different regime. Beyond fixed-boundary IVPs, the package also solves **endogenous free-boundary problems** — including problems requiring several independently-integrated legs glued together by value-matching / smooth-pasting conditions, such as real-options entry/exit models — entirely through its own `solve` interface.

The package provides:

- Initial value problem (IVP) switching solvers
- Event-driven regime transitions
- Pre-switch and post-switch trajectory handling
- ContinuousCallback-based switching
- Endogenous free-boundary solvers (value-matching / smooth-pasting)
- Multi-leg free-boundary solvers for problems with several independently-integrated regimes and thresholds
- Compatibility with the SciML ecosystem
- Infrastructure for shooting-based boundary value methods

---

# Features

- Endogenous switching using user-defined conditions
- Separate dynamics before and after switching
- Continuous event detection
- Free-boundary (value-matching / smooth-pasting) solving via `mode = :free_boundary`
- Multi-leg free-boundary solving via `Leg` and `mode = :multi_leg`, for problems glued together across several separately-integrated regimes (e.g. entry/exit, capacity choice)
- Type-stable solution containers
- DifferentialEquations.jl integration
- The full nonlinear solve — including algorithm selection — stays inside the package; no need to depend on `NonlinearSolve.jl` directly in your own code
- Simple and extensible API
- Built on the SciML ecosystem

---

# Installation

```julia
using Pkg
Pkg.add("SwitchingControl")
```

`SwitchingControl.jl` is available through the Julia General Registry.

---

# Quick Start

## Basic Switching Example

```julia
using SwitchingControl

# ============================================================
# PARAMETERS
# ============================================================

p = nothing

# ============================================================
# PRE-SWITCH DYNAMICS
# ============================================================

function dynamics_before!(du, u, p, t)

    x = u[1]
    A = u[2]

    du[1] = 0.2 * x
    du[2] = 1.0

    return nothing
end

# ============================================================
# POST-SWITCH DYNAMICS
# ============================================================

function dynamics_after!(du, u, p, t)

    x = u[1]
    A = u[2]

    du[1] = -0.3 * x
    du[2] = -0.5

    return nothing
end

# ============================================================
# SWITCHING CONDITION
# ============================================================

function switching_condition(u, p, t)

    A = u[2]

    return A - 5.0
end

# ============================================================
# BOUNDARY CONDITIONS
# ============================================================

bc = BoundaryConditions(
    [1.0, 0.0],
    nothing
)

# ============================================================
# PROBLEM
# ============================================================

prob = SwitchingProblem(
    dynamics_before!,
    dynamics_after!,
    switching_condition,
    bc,
    p,
    (0.0, 50.0),
    1
)

# ============================================================
# SOLVE
# ============================================================

sol = SwitchingControl.solve(prob)

println("Switched?      : ", sol.sol.switched)
println("Switch time    : ", sol.sol.tau)
println("Converged?     : ", sol.sol.converged)
```

---

# Solution Structure

The solver returns a `SwitchingSolutionView` object.

Main fields:

```julia
sol.sol.tau
```

Switching time.

```julia
sol.sol.switched
```

Whether switching occurred.

```julia
sol.sol.solution_before
```

ODE solution before switching.

```julia
sol.sol.solution_after
```

ODE solution after switching.

```julia
sol.sol.converged
```

Solver convergence status.

---

# Example: Accessing Trajectories

## Pre-switch states

```julia
u_pre = sol.sol.solution_before.u
```

## Post-switch states

```julia
u_post = sol.sol.solution_after.u
```

## Time grids

```julia
t_pre = sol.sol.solution_before.t
t_post = sol.sol.solution_after.t
```

---

# Switching Logic

The package uses `ContinuousCallback` from DifferentialEquations.jl.

A switch occurs when:

```math
g(u, p, t) = 0
```

where `g` is the user-defined switching condition.

The solver:

1. Integrates the pre-switch system
2. Detects the switching event
3. Terminates integration
4. Restarts using post-switch dynamics
5. Returns the combined switching solution

---

# Free-Boundary Problems (Value-Matching / Smooth-Pasting)

Many optimal-stopping and real-options problems don't have a fixed switching time — the switching *threshold itself* is unknown and pinned down by economic conditions such as value-matching and smooth-pasting. `mode = :free_boundary` solves exactly this: given a function that builds the initial state from free parameters `θ` and a residual evaluated on the resulting switching solution, the package handles the entire outer nonlinear solve internally.

Example — McDonald & Siegel (1986) irreversible investment:

```julia
using SwitchingControl

r, δ, σ, I = 0.05, 0.04, 0.20, 1.0
p = (r = r, δ = δ, σ = σ, I = I)

# Continuation region
function dynamics_before!(du, u, p, V)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*V*Fp) / (p.σ^2 * V^2)
    return nothing
end

# Exercised region: F(V) = V - I
function dynamics_after!(du, u, p, V)
    du[1] = u[2]
    du[2] = 0.0
    return nothing
end

switching_condition(u, p, V) = u[2] - 1.0   # smooth pasting: F'(V) = 1

const ε = 1e-3
make_u0(θ) = [θ[1]*ε^0.5, θ[1]*0.5*ε^(0.5 - 1)]  # placeholder exponent, see examples/mcdonald_siegel.jl

function residual!(F, sol, θ, p)
    Vswitch = sol.sol.tau
    Fswitch = sol.sol.solution_before.u[end][1]
    F[1] = Fswitch - (Vswitch - p.I)   # value matching
    return nothing
end

bc = BoundaryConditions([0.0, 0.0], nothing)  # placeholder; make_u0 supplies the real u0

prob = SwitchingProblem(dynamics_before!, dynamics_after!, switching_condition, bc, p, (ε, 10.0*I), 1)

nlsol = SwitchingControl.solve(
    prob;
    mode = :free_boundary,
    guess = [0.25],
    make_u0 = make_u0,
    residual! = residual!
)
```

See `examples/mcdonald_siegel.jl` for the full, runnable version (with the correct closed-form exponent) plus validation against the known closed-form threshold.

---

# Multi-Leg Free-Boundary Problems

Some free-boundary problems — most notably entry/exit and capacity-choice models — can't be reduced to a single IVP with one switch. They require **several independently-integrated legs** (e.g. an idle-regime ODE and an active-regime ODE, sometimes integrated in opposite directions over different ranges) glued together by a *joint* residual over shared free parameters. `mode = :multi_leg` (or the `solve(legs::Vector{Leg}, guess; ...)` convenience form) handles this natively:

- Each `Leg` specifies its own `dynamics!`, `tspan`, `parameters`, `make_u0(θ)`, and `switching_condition(θ)`.
- On every nonlinear-solver iterate, every leg is re-solved for the current `θ`, and the resulting solutions are handed to your `residual!(F, sols, θ, p)` in the same order as `legs`.
- If any leg's guessed threshold falls outside its integrated range, the package fills a penalty residual instead of throwing, so globalized algorithms can backtrack rather than crash.
- The nonlinear algorithm and autodiff scheme can be selected with plain symbols (`:trust_region`, `:newton`, `:levenberg_marquardt`; `:central`, `:forward_diff`, `:forward`) — no need to `using NonlinearSolve` yourself.

Example — Dixit (1989) entry/exit under uncertainty:

```julia
using SwitchingControl

# ... dynamics_v0!, dynamics_v1! defined as in examples/dixit.jl ...

v0_leg = Leg(
    dynamics! = dynamics_v0!,
    tspan = (ε, P_upper),
    parameters = p,
    make_u0 = θ -> [θ[1]*ε^β1, θ[1]*β1*ε^(β1 - 1)],
    switching_condition = θ -> ((u, p, P) -> P - θ[3])
)

v1_leg = Leg(
    dynamics! = dynamics_v1!,
    tspan = (P_far, P_lower),   # integrated backward
    parameters = p,
    make_u0 = θ -> [θ[2]*P_far^β2 + P_far/p.δ - p.w/p.r,
                    θ[2]*β2*P_far^(β2 - 1) + 1/p.δ],
    switching_condition = θ -> ((u, p, P) -> P - θ[4])
)

function dixit_residual!(F, sols, θ, p)
    sol_v0, sol_v1 = sols[1], sols[2]
    V0_at_PH = sol_v0.sol.solution_before.u[end]
    V1_at_PL = sol_v1.sol.solution_before.u[end]
    V0_at_PL = sol_v0.sol.solution_before(θ[4])
    V1_at_PH = sol_v1.sol.solution_before(θ[3])

    F[1] = V1_at_PH[1] - p.I - V0_at_PH[1]   # entry value-matching
    F[2] = V1_at_PH[2] - V0_at_PH[2]         # entry smooth-pasting
    F[3] = V0_at_PL[1] - (V1_at_PL[1] - p.E) # exit value-matching
    F[4] = V0_at_PL[2] - V1_at_PL[2]         # exit smooth-pasting
    return nothing
end

nlsol = SwitchingControl.solve(
    [v0_leg, v1_leg], [32.0, 2.8, 1.2, 0.85];
    residual! = dixit_residual!,
    autodiff = :central,
    ivp_kwargs = (abstol = 1e-12, reltol = 1e-12)
)
```

See `examples/dixit.jl` for the full, runnable version, validated against the closed-form entry/exit thresholds.

---

# Shooting Solver

The package also contains infrastructure for shooting-based boundary value methods.

Example:

```julia
result = SwitchingControl.solve(
    prob;
    mode = :shooting,
    guess = [0.2, 0.1]
)
```

Current shooting support is experimental and works best for smooth, well-conditioned systems.

---

# Choosing a Solve Mode

| Mode              | Use when...                                                                 |
|-------------------|------------------------------------------------------------------------------|
| `:ivp` (default)  | You know the switching *condition*, and just want to integrate through it   |
| `:shooting`        | Classic two-point BVP: match a fixed terminal state at a fixed final time   |
| `:free_boundary`   | One IVP whose switching *threshold* is itself unknown (value-matching / smooth-pasting) |
| `:multi_leg`       | Several independently-integrated legs (e.g. idle/active regimes) glued by a joint residual over shared free parameters |

---

# Limitations

Current limitations of the package:

- Each individual `SwitchingProblem`/`Leg` supports only a single endogenous switching event per solve
- Shooting-based methods are currently experimental
- `:multi_leg` covers problems expressible as several separately-integrated legs sharing free parameters; general arbitrary-length regime *sequencing* within one continuous trajectory (e.g. idle → active → idle → active, endogenously repeating) is not yet implemented
- Primarily designed for IVP and free-boundary switching problems
- PDE and HJB functionality are planned but not yet available

---

# Package Structure

```text
SwitchingControl.jl/
│
├── src/
│   ├── SwitchingControl.jl
│   ├── Types.jl
│   ├── Solution.jl
│   └── solver.jl
│
├── test/
│   └── runtests.jl
│
├── examples/
│   ├── basic_switching.jl
│   ├── mcdonald_siegel.jl
│   └── dixit.jl
│
├── Project.toml
└── README.md
```

---

# Running Tests

From the Julia REPL:

```julia
using Pkg

Pkg.test("SwitchingControl")
```

---

# Running the Examples

From the Julia REPL:

```julia
include("examples/basic_switching.jl")
include("examples/mcdonald_siegel.jl")
include("examples/dixit.jl")
```

---

# SciML Ecosystem

`SwitchingControl.jl` is built on top of:

- DifferentialEquations.jl
- SciMLBase.jl
- NonlinearSolve.jl
- BoundaryValueDiffEq.jl

`NonlinearSolve.jl` and its AD backends are used internally to power `:free_boundary` and `:multi_leg` mode — you don't need to depend on them directly; `autodiff`/`nl_alg` accept plain symbols instead.

---

# Planned Features

- Multiple switching events within a single continuous trajectory
- Hybrid system support
- Additional optimal control interfaces
- Economic control templates
- PDE-based extensions
- Visualization utilities
- Improved shooting methods
- Automatic continuity enforcement
- General regime sequence support

---

# Applications

Potential use cases include:

- Optimal stopping problems
- Entry/exit and capacity-choice models under uncertainty
- Retirement models
- Pandemic economics
- Hybrid dynamical systems
- Mechanical switching systems
- Threshold control problems
- State-triggered policy rules
- Dynamic optimization

---

# License

MIT License
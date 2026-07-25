# SwitchingControl.jl

![CI](https://github.com/vishalkumarrai-1994/SwitchingControl.jl/actions/workflows/ci.yml/badge.svg)


A Julia framework for solving dynamic optimization and switching-control problems through event-driven numerical integration.

`SwitchingControl.jl` is designed for optimal control, economics, engineering, and hybrid dynamical systems where trajectories evolve under one regime until a switching condition is triggered, after which the system transitions to a different regime.

The package provides:

- Initial value problem (IVP) switching solvers
- Event-driven regime transitions
- Pre-switch and post-switch trajectory handling
- ContinuousCallback-based switching
- Compatibility with the SciML ecosystem
- Infrastructure for shooting-based boundary value methods

---

# Features

- Endogenous switching using user-defined conditions
- Separate dynamics before and after switching
- Continuous event detection
- Type-stable solution containers
- DifferentialEquations.jl integration
- Simple and extensible API
- Built on the SciML ecosystem

---

# Installation

# Installation

```julia
using Pkg
Pkg.add("SwitchingControl")
```

`SwitchingControl.jl` is available through the Julia General Registry.
----

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

# Limitations

Current limitations of the package:

- Only a single endogenous switching event is supported per solve
- Shooting-based methods are currently experimental
- Multiple-regime sequencing is not yet implemented
- Primarily designed for IVP switching problems
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
│   └── basic_switching.jl
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

# Running the Example

From the Julia REPL:

```julia
include("examples/basic_switching.jl")
```

---

# SciML Ecosystem

`SwitchingControl.jl` is built on top of:

- DifferentialEquations.jl
- SciMLBase.jl
- NonlinearSolve.jl
- BoundaryValueDiffEq.jl

---

# Planned Features

- Multiple switching events
- Hybrid system support
- Additional optimal control interfaces
- Economic control templates
- PDE-based extensions
- Visualization utilities
- Improved shooting methods
- Automatic continuity enforcement
- Regime sequence support

---

# Applications

Potential use cases include:

- Optimal stopping problems
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

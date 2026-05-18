module Types

export SwitchingProblem, SwitchingSolution, BoundaryConditions

# -----------------------------
# Boundary conditions
# -----------------------------
"""
Stores initial or boundary conditions for switching problems.

Fields typically contain initial values and optional terminal conditions.
"""
struct BoundaryConditions{U, V}
    initial::U
    terminal::V
end


# -----------------------------
# Switching problem
# -----------------------------
"""
Defines a dynamical system with endogenous regime switching.

The system evolves under `f_before!` until the switching condition
becomes zero, after which dynamics transition to `f_after!`.
"""
struct SwitchingProblem{F1, F2, S, BC, P, TS}
    dynamics_before!::F1
    dynamics_after!::F2

    switching_condition::S

    boundary_conditions::BC
    parameters::P
    tspan::TS

    max_switches::Int
end


# -----------------------------
# Raw solution container (internal)
# -----------------------------
"""
Stores the full solution and metadata for a switching problem.
"""
struct SwitchingSolution{TT, SB, SA}
    tau::TT
    solution_before::SB
    solution_after::SA
    switched::Bool
    converged::Bool
    residual_norm::Float64
end

end
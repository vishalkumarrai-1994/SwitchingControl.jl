module Types

export SwitchingProblem, SwitchingSolution, BoundaryConditions

# -----------------------------
# Boundary conditions
# -----------------------------
struct BoundaryConditions{U, V}
    initial::U
    terminal::V
end


# -----------------------------
# Switching problem
# -----------------------------
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
struct SwitchingSolution{TT, SB, SA}
    tau::TT
    solution_before::SB
    solution_after::SA
    switched::Bool
    converged::Bool
    residual_norm::Float64
end

end
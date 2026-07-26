module Types

export AbstractSwitchingProblem, SwitchingProblem, SwitchingSolution, BoundaryConditions

abstract type AbstractSwitchingProblem end

struct BoundaryConditions{U, V}
    initial::U
    terminal::V
end

struct SwitchingProblem{F1, F2, S, BC, P, TS} <: AbstractSwitchingProblem
    dynamics_before!::F1
    dynamics_after!::F2
    switching_condition::S
    boundary_conditions::BC
    parameters::P
    tspan::TS
    max_switches::Int
end

struct SwitchingSolution{TT, SB, SA}
    tau::TT
    solution_before::SB
    solution_after::SA
    switched::Bool
    converged::Bool
    residual_norm::Float64
end

end
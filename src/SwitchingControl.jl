module SwitchingControl

using DifferentialEquations
using SciMLBase
using NonlinearSolve

include("types.jl")
include("solution.jl")
include("solver.jl")

using .Types
using .Solution

export SwitchingProblem,
       SwitchingSolution,
       SwitchingSolutionView,
       BoundaryConditions,
       solve

end
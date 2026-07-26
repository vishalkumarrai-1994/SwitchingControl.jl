module SwitchingControl

using DifferentialEquations, SciMLBase, NonlinearSolve, Polynomials, LinearAlgebra

function solve end   # generic function — both solver files add methods

include("types.jl");        using .Types
include("MarkovTypes.jl");  using .MarkovTypes
include("solution.jl");     using .Solution
include("solver.jl")
include("MarkovSolver.jl")

export SwitchingProblem, SwitchingSolution, SwitchingSolutionView,
       MarkovSwitchingProblem, MarkovSwitchingSolution, MarkovSwitchingSolutionView,
       GBMRegime, GenericRegime, PutPayoff, GenericPayoff,
       BoundaryConditions, AbstractSwitchingProblem, solve
export Leg

end
module MarkovTypes

export AbstractRegimeDynamics, GBMRegime, GenericRegime,
       AbstractPayoff, PutPayoff, GenericPayoff,
       MarkovSwitchingProblem, MarkovSwitchingSolution

using ..Types: AbstractSwitchingProblem

# ------------------------------------------------------------
# Regime dynamics traits
# ------------------------------------------------------------
abstract type AbstractRegimeDynamics end

"""
Constant-coefficient GBM regime: dX = μX dt + σX dW.
Enables the exact closed-form solver (Guo & Zhang 2004) for S=2.
"""
struct GBMRegime <: AbstractRegimeDynamics
    μ::Float64
    σ::Float64
end

"""
General regime dynamics via user-supplied (x,p) -> (drift, diffusion).
Routed to the numeric coupled-BVP solver.
"""
struct GenericRegime{F} <: AbstractRegimeDynamics
    f::F
end

# ------------------------------------------------------------
# Payoff traits
# ------------------------------------------------------------
abstract type AbstractPayoff end

struct PutPayoff <: AbstractPayoff
    K::Float64
end
(p::PutPayoff)(x) = max(p.K - x, 0.0)

struct GenericPayoff{F} <: AbstractPayoff
    f::F
end
(p::GenericPayoff)(x) = p.f(x)

# ------------------------------------------------------------
# Problem / solution containers
# ------------------------------------------------------------
"""
N regimes coupled through a Markov generator, each with its own state
dynamics and an endogenous free boundary where optimal stopping occurs.
Unlike `SwitchingProblem`, regimes are not temporally sequenced — all N
value functions are solved jointly over the same state domain.
"""
struct MarkovSwitchingProblem{D<:AbstractRegimeDynamics, QT, PT<:AbstractPayoff, PM, DM} <: AbstractSwitchingProblem
    regime_dynamics::Vector{D}
    generator::QT
    payoff::Vector{PT}
    discount::Float64
    domain::DM
    parameters::PM
    n_regimes::Int
end

struct MarkovSwitchingSolution{B, V}
    boundaries::B
    value_functions::V
    converged::Bool
    residual_norm::Float64
end

end
module MarkovSolver

using ..MarkovTypes
using ..Solution: MarkovSwitchingSolutionView
using LinearAlgebra, Polynomials, NonlinearSolve, SciMLBase

import ..SwitchingControl: solve

# ============================================================
# Exact algebraic path: S=2, GBM regimes, put payoff
# (Guo & Zhang 2004, "Closed-Form Solutions for Perpetual
# American Put Options with Regime Switching", eqs. 8-14)
# ============================================================
function solve(prob::MarkovSwitchingProblem{GBMRegime,QT,PutPayoff};
               guess = nothing, kwargs...) where {QT}

    prob.n_regimes == 2 || error(
        "Exact closed-form solver applies only for S=2 regimes " *
        "(Guo & Zhang 2004, Remark 3.5: closed form exists iff S=2). " *
        "Use GenericRegime dynamics with n_regimes > 2 for the numeric solver."
    )

    r  = prob.discount
    K  = prob.payoff[1].K
    μ  = Dict(1 => prob.regime_dynamics[1].μ, 2 => prob.regime_dynamics[2].μ)
    σ  = Dict(1 => prob.regime_dynamics[1].σ, 2 => prob.regime_dynamics[2].σ)
    λ  = Dict(1 => prob.generator[1,2],       2 => prob.generator[2,1])

    guess === nothing && (guess = [0.4K, 0.6K])

    for (active, pinned) in ((1,2), (2,1))
        result = _try_solve_case(r, μ, σ, λ, K, active, pinned, guess)
        result !== nothing && return _to_view(result, prob)
    end
    error("No consistent threshold ordering found for the given parameters.")
end

function _characteristic_roots(r, μ, σ, λ)
    g1c = [λ[1]+r, -(μ[1]-0.5σ[1]^2), -0.5σ[1]^2]
    g2c = [λ[2]+r, -(μ[2]-0.5σ[2]^2), -0.5σ[2]^2]
    prodpoly = Polynomial(g1c) * Polynomial(g2c) - λ[1]*λ[2]
    rts = roots(prodpoly)
    real_rts = sort(real.(filter(z -> abs(imag(z)) < 1e-8, rts)))
    length(real_rts) >= 2 || error("Expected >=2 real characteristic roots; got $(real_rts)")
    return real_rts[1], real_rts[2]   # β1 < β2 < 0
end

_g1(β, r, μ, σ, λ) = λ[1] + r - (μ[1]-0.5σ[1]^2)*β - 0.5σ[1]^2*β^2

function _try_solve_case(r, μ, σ, λ, K, active, pinned, guess)
    β1, β2 = _characteristic_roots(r, μ, σ, λ)
    l1 = _g1(β1, r, μ, σ, λ) / λ[1]
    l2 = _g1(β2, r, μ, σ, λ) / λ[1]

    a  = 0.5σ[active]^2
    b_ = μ[active] - 0.5σ[active]^2
    c_ = -(r + λ[active])
    disc = b_^2 - 4a*c_
    disc < 0 && return nothing
    γ1 = (-b_ + sqrt(disc)) / (2a)
    γ2 = (-b_ - sqrt(disc)) / (2a)
    c0 = λ[active]*K / (r+λ[active])
    c1 = -λ[active] / (r+λ[active]-μ[active])
    φ(x)  = c0 + c1*x
    φp(x) = c1

    function solve_A(x_hi)
        M = pinned == 1 ?
            [x_hi^β1        x_hi^β2;
             β1*x_hi^β1     β2*x_hi^β2] :
            [l1*x_hi^β1        l2*x_hi^β2;
             β1*l1*x_hi^β1     β2*l2*x_hi^β2]
        M \ [K - x_hi, -x_hi]
    end

    function V_active_outer(x, A1, A2)
        if active == 1
            val = A1*x^β1 + A2*x^β2
            der = β1*A1*x^β1 + β2*A2*x^β2
        else
            val = l1*A1*x^β1 + l2*A2*x^β2
            der = β1*l1*A1*x^β1 + β2*l2*A2*x^β2
        end
        return val, der
    end

    function solve_C(x_lo)
        M = [x_lo^γ1        x_lo^γ2;
             γ1*x_lo^γ1     γ2*x_lo^γ2]
        M \ [K - x_lo - φ(x_lo), -x_lo - x_lo*φp(x_lo)]
    end

    function residual!(F, xs, p)
        x_lo, x_hi = xs
        if x_lo <= 0 || x_hi <= 0 || x_lo >= x_hi
            F .= 1e6; return
        end
        A1, A2 = solve_A(x_hi)
        C1, C2 = solve_C(x_lo)
        val_a, der_a = V_active_outer(x_hi, A1, A2)
        val_m = C1*x_hi^γ1 + C2*x_hi^γ2 + φ(x_hi)
        der_m = γ1*C1*x_hi^γ1 + γ2*C2*x_hi^γ2 + x_hi*φp(x_hi)
        F[1] = val_a - val_m
        F[2] = der_a - der_m
    end

    nlprob = NonlinearProblem(residual!, guess, nothing)
    nlsol  = NonlinearSolve.solve(nlprob, NewtonRaphson())
    !SciMLBase.successful_retcode(nlsol) && return nothing

    x_lo, x_hi = nlsol.u
    (x_lo <= 0 || x_hi <= 0 || x_lo >= x_hi) && return nothing

    boundaries = active == 1 ? [x_lo, x_hi] : [x_hi, x_lo]  # map back to regime 1,2 order
    A1, A2 = solve_A(x_hi)
    C1, C2 = solve_C(x_lo)

    (; boundaries, A1, A2, C1, C2, β1, β2, γ1, γ2, l1, l2, φ, active, pinned,
       resid = norm(nlsol.resid))
end

function _to_view(result, prob::MarkovSwitchingProblem)
    raw = MarkovSwitchingSolution(
        result.boundaries,
        result,  # keep full algebraic pieces for downstream evaluation
        true,
        result.resid
    )
    return MarkovSwitchingSolutionView(raw)
end

# ============================================================
# Generic fallback: N regimes, arbitrary dynamics/payoff
# (numeric coupled-BVP path — TODO, next milestone)
# ============================================================
function solve(prob::MarkovSwitchingProblem; kwargs...)
    error("Generic N-regime numeric solver not yet implemented — " *
          "currently only GBMRegime + PutPayoff with n_regimes==2 is supported.")
end

end
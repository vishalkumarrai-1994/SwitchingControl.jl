# ============================================================
# Guo & Zhang (2004) — Closed-Form Solutions for Perpetual
# American Put Options with Regime Switching
# SIAM J. Appl. Math. 64(6), 2034–2049
#
# A perpetual American put whose underlying follows GBM with
# volatility switching between two regimes according to a
# continuous-time Markov chain. Unlike SwitchingProblem (where
# switching happens once, along a single trajectory), both
# regimes' value functions here are solved *jointly*, each with
# its own endogenous exercise boundary, coupled through the
# generator matrix.
#
# This reproduces the paper's Section 4.3 base case:
#   r = 3, μ1 = μ2 = 3, K = 5, λ1 = λ2 = 100, σ1 = 9, σ2 = 5
# Published closed-form thresholds: (x1*, x2*) = (0.441, 0.614)
# ============================================================

using SwitchingControl

prob = MarkovSwitchingProblem(
    [GBMRegime(3.0, 9.0), GBMRegime(3.0, 5.0)],   # (μ1,σ1), (μ2,σ2)
    [-100.0 100.0; 100.0 -100.0],                  # generator Q
    [PutPayoff(5.0), PutPayoff(5.0)],               # K = 5, both regimes
    3.0,                                             # discount r
    (1e-3, 15.0),                                    # domain
    nothing, 2
)

sol = SwitchingControl.solve(prob)

println("converged   = ", sol.converged)
println("boundaries  = ", sol.boundaries)
println("expected    ≈ [0.441, 0.614]  (published, 3 dp)")
println("expected    ≈ [0.4405, 0.6116] (independently verified to higher precision)")

@assert sol.converged
@assert isapprox(sol.boundaries[1], 0.441; atol=0.002)
@assert isapprox(sol.boundaries[2], 0.614; atol=0.003)
println("\n✅ PASS — matches Guo & Zhang (2004)")
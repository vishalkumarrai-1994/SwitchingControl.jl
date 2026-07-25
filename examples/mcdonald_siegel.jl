using SwitchingControl

# ------------------------------------------------------------
# Parameters (McDonald & Siegel, 1986)
# ------------------------------------------------------------
r = 0.05   # discount rate
δ = 0.04   # payout / cost-of-delay yield on V
σ = 0.20   # volatility
I = 1.0    # investment cost

# Closed-form beta and threshold — for validation only.
β = 0.5 - (r - δ)/σ^2 + sqrt(((r - δ)/σ^2 - 0.5)^2 + 2r/σ^2)
Vstar_closed_form = β/(β - 1) * I

p = (r = r, δ = δ, σ = σ, I = I)

# ------------------------------------------------------------
# Treat V as the package's independent variable "t".
# State u = [F(V), F'(V)].
# ------------------------------------------------------------

# Continuation region: 0.5σ²V²F'' + (r-δ)VF' - rF = 0
function dynamics_before!(du, u, p, V)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*V*Fp) / (p.σ^2 * V^2)
    return nothing
end

# Exercised region: F(V) = V - I  =>  F' = 1, F'' = 0
function dynamics_after!(du, u, p, V)
    du[1] = u[2]
    du[2] = 0.0
    return nothing
end

# Smooth pasting: switch exactly when F'(V) = 1
switching_condition(u, p, V) = u[2] - 1.0

const ε = 1e-3

# ------------------------------------------------------------
# Free-boundary mode: package now owns the outer solve.
#   θ = [A], the local-expansion coefficient F(V) ≈ A·V^β near V=0
#   residual = value-matching, evaluated AT THE SWITCH POINT
# ------------------------------------------------------------

make_u0(θ) = [θ[1]*ε^β, θ[1]*β*ε^(β - 1)]

function residual!(F, sol, θ, p)
    @assert sol.sol.switched "smooth pasting (F'=1) never reached — widen tspan"
    Vswitch = sol.sol.tau
    Fswitch = sol.sol.solution_before.u[end][1]
    F[1] = Fswitch - (Vswitch - p.I)   # value matching
    return nothing
end

bc = BoundaryConditions([0.0, 0.0], nothing)  # placeholder; make_u0 supplies the real u0

prob = SwitchingProblem(
    dynamics_before!,
    dynamics_after!,
    switching_condition,
    bc, p,
    (ε, 10.0*I),
    1
)

nlsol = SwitchingControl.solve(
    prob;
    mode = :free_boundary,
    guess = [0.25],
    make_u0 = make_u0,
    residual! = residual!
)

A_star = nlsol.u[1]

# One more solve at the converged A to report V* directly.
u0_star = make_u0([A_star])
final_prob = SwitchingProblem(dynamics_before!, dynamics_after!, switching_condition,
                               BoundaryConditions(u0_star, nothing), p, (ε, 10.0*I), 1)
final_sol = SwitchingControl.solve(final_prob)
Vstar_numeric = final_sol.sol.tau

println("Closed-form V* : ", Vstar_closed_form)
println("Numerical  V*  : ", Vstar_numeric)
println("Relative error : ", abs(Vstar_numeric - Vstar_closed_form)/Vstar_closed_form)


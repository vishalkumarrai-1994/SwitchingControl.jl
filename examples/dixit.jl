using SwitchingControl

# ------------------------------------------------------------
# Parameters (Dixit, 1989, "Entry and Exit Decisions Under
# Uncertainty", JPE)
# ------------------------------------------------------------
r = 0.05   # discount rate
δ = 0.02   # payout / cost-of-delay yield on P
σ = 0.20   # volatility
w = 1.0    # fixed operating cost while active
I = 0.7    # entry (sunk) cost
E = 0.3    # exit (sunk) cost

# Closed-form thresholds — for validation only.
A2 = 0.5*σ^2
A1 = (r - δ) - 0.5*σ^2
A0 = -r
disc = sqrt(A1^2 - 4*A2*A0)
β1 = (-A1 + disc) / (2*A2)   # > 1
β2 = (-A1 - disc) / (2*A2)   # < 0

P_H_closed_form = 1.2701531078
P_L_closed_form = 0.7987448046

p = (r = r, δ = δ, σ = σ, w = w, I = I, E = E)

# ------------------------------------------------------------
# Treat P as the package's independent variable "t".
# State u = [F(P), F'(P)], one leg per regime:
#   V0 (idle):   0.5σ²P²F'' + (r-δ)PF' - rF = 0
#   V1 (active): 0.5σ²P²F'' + (r-δ)PF' - rF = -(P - w)
# ------------------------------------------------------------

function dynamics_v0!(du, u, p, P)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*P*Fp) / (p.σ^2 * P^2)
    return nothing
end

function dynamics_v1!(du, u, p, P)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*P*Fp - 2*(P - p.w)) / (p.σ^2 * P^2)
    return nothing
end

const ε       = 0.10   # start of the V0 leg -- moved away from the P=0
                        # singularity; a*P^β1 is exact for all P>0, so this
                        # loses nothing as long as it stays below P_L.
const P_upper = 5.0    # end of the V0 leg's domain (must exceed the true P_H)
const P_far   = 50.0   # start of the V1 leg (large P, asymptotic regime)
const P_lower = 0.01   # end of the V1 leg's domain (must be below the true P_L)

# ------------------------------------------------------------
# Multi-leg free-boundary mode: package owns the outer solve.
#   θ = (a, b, P_H, P_L)
#     a, b   -- local-expansion / asymptotic coefficients
#     P_H    -- entry threshold (V0 leg's target)
#     P_L    -- exit threshold  (V1 leg's target)
#   residual = entry value-matching, entry smooth-pasting,
#              exit value-matching,  exit smooth-pasting
# ------------------------------------------------------------

v0_leg = Leg(
    dynamics! = dynamics_v0!,
    tspan = (ε, P_upper),
    parameters = p,
    make_u0 = θ -> begin
        a = θ[1]
        [a*ε^β1, a*β1*ε^(β1 - 1)]
    end,
    switching_condition = θ -> begin
        P_H_guess = θ[3]
        (u, p, P) -> P - P_H_guess
    end
)

v1_leg = Leg(
    dynamics! = dynamics_v1!,
    tspan = (P_far, P_lower),      # backward
    parameters = p,
    make_u0 = θ -> begin
        b = θ[2]
        [b*P_far^β2 + P_far/p.δ - p.w/p.r,
         b*β2*P_far^(β2 - 1) + 1/p.δ]
    end,
    switching_condition = θ -> begin
        P_L_guess = θ[4]
        (u, p, P) -> P - P_L_guess
    end
)

function dixit_residual!(F, sols, θ, p)
    a, b, P_H, P_L = θ

    sol_v0, sol_v1 = sols[1], sols[2]

    V0_at_PH = sol_v0.sol.solution_before.u[end]     # exact terminal state
    V1_at_PL = sol_v1.sol.solution_before.u[end]     # exact terminal state

    V0_at_PL = sol_v0.sol.solution_before(P_L)       # interpolated: P_L ∈ (ε, P_H)
    V1_at_PH = sol_v1.sol.solution_before(P_H)       # interpolated: P_H ∈ (P_L, P_far)

    F[1] = V1_at_PH[1] - p.I - V0_at_PH[1]           # entry value-matching
    F[2] = V1_at_PH[2] - V0_at_PH[2]                 # entry smooth-pasting
    F[3] = V0_at_PL[1] - (V1_at_PL[1] - p.E)         # exit value-matching
    F[4] = V0_at_PL[2] - V1_at_PL[2]                 # exit smooth-pasting

    return nothing
end

guess = [32.0, 2.8, 1.2, 0.85]

nlsol = SwitchingControl.solve(
    [v0_leg, v1_leg], guess;
    residual! = dixit_residual!,
    autodiff = :central,
    ivp_kwargs = (abstol = 1e-12, reltol = 1e-12),
    nlsolve_kwargs = (abstol = 1e-9, maxiters = 100)
)

a_star, b_star, P_H_star, P_L_star = nlsol.u

println("Closed-form P_H : ", P_H_closed_form)
println("Numerical  P_H  : ", P_H_star)
println("Relative error  : ", abs(P_H_star - P_H_closed_form)/P_H_closed_form)
println()
println("Closed-form P_L : ", P_L_closed_form)
println("Numerical  P_L  : ", P_L_star)
println("Relative error  : ", abs(P_L_star - P_L_closed_form)/P_L_closed_form)
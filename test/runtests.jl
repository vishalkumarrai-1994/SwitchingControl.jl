using Test
using LinearAlgebra
using SciMLBase
using SwitchingControl

# ============================================================
# SIMPLE TEST PROBLEM
# ============================================================

function f_before!(dy, y, p, t)
    x, A = y
    dy[1] = 1.0
    dy[2] = -0.5 * A
    return nothing
end

function f_after!(dy, y, p, t)
    x, A = y
    dy[1] = -0.25
    dy[2] = 0.1 * A
    return nothing
end

function switching_condition(y, p, t)
    return y[1] - 1.0
end

bc = BoundaryConditions([0.0, 10.0], nothing)

prob = SwitchingProblem(
    f_before!,
    f_after!,
    switching_condition,
    bc,
    nothing,
    (0.0, 5.0),
    1
)

@testset "SwitchingControl IVP Tests" begin

    sol = SwitchingControl.solve(prob; mode = :ivp)

    @test sol.sol.converged == true
    @test sol.sol.switched == true
    @test !isnothing(sol.sol.tau)

    @test sol.sol.solution_before !== nothing
    @test sol.sol.solution_after !== nothing

    dims_pre = unique(length(u) for u in sol.sol.solution_before.u)
    dims_post = unique(length(u) for u in sol.sol.solution_after.u)

    @test dims_pre == [2]
    @test dims_post == [2]

    u_before = sol.sol.solution_before.u[end]
    u_after  = sol.sol.solution_after.u[1]

    @test norm(u_after .- u_before) < 1e-12
    @test abs(u_before[1] - 1.0) < 1e-6

    @test all(isfinite, sol.sol.solution_after.u[end])

    @test length(sol.sol.solution_before.u) > 1
    @test length(sol.sol.solution_after.u) > 1
end


# ============================================================
# McDONALD-SIEGEL (1986) FREE-BOUNDARY / SMOOTH-PASTING PROBLEM
#
# Helper functions are prefixed `ms_` so they don't collide with
# (or silently overwrite) same-named functions in other testsets.
# ============================================================

const MS_r = 0.05     # discount rate
const MS_δ = 0.04     # payout / cost-of-delay yield on V
const MS_σ = 0.20     # volatility
const MS_I0 = 1.0     # investment cost (renamed from `I` — LinearAlgebra exports `I`)
const MS_ε = 1e-3     # start just above the coordinate singularity at V=0

# Closed-form beta and threshold, used only to validate the numerics.
const MS_β = 0.5 - (MS_r - MS_δ)/MS_σ^2 +
             sqrt(((MS_r - MS_δ)/MS_σ^2 - 0.5)^2 + 2*MS_r/MS_σ^2)
const MS_Vstar_closed_form = MS_β/(MS_β - 1) * MS_I0

const MS_p = (r = MS_r, δ = MS_δ, σ = MS_σ, I = MS_I0)

# Continuation region: 0.5σ²V²F'' + (r-δ)VF' - rF = 0
function ms_dynamics_before!(du, u, p, V)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*V*Fp) / (p.σ^2 * V^2)
    return nothing
end

# Exercised region: F(V) = V - I  =>  F' = 1, F'' = 0
function ms_dynamics_after!(du, u, p, V)
    du[1] = u[2]
    du[2] = 0.0
    return nothing
end

# Smooth pasting: switch exactly when F'(V) = 1
ms_switching_condition(u, p, V) = u[2] - 1.0

# θ = [A], the local-expansion coefficient F(V) ≈ A·V^β near V=0
ms_make_u0(θ) = [θ[1]*MS_ε^MS_β, θ[1]*MS_β*MS_ε^(MS_β - 1)]

# Value-matching residual, evaluated AT THE SWITCH POINT.
function ms_residual!(F, sol, θ, p)
    @assert sol.sol.switched "smooth pasting (F'=1) never reached — widen tspan"
    Vswitch = sol.sol.tau
    Fswitch = sol.sol.solution_before.u[end][1]
    F[1] = Fswitch - (Vswitch - p.I)
    return nothing
end

@testset "McDonald-Siegel free boundary" begin

    ms_bc = BoundaryConditions([0.0, 0.0], nothing)  # placeholder; ms_make_u0 supplies real u0

    ms_prob = SwitchingProblem(
        ms_dynamics_before!,
        ms_dynamics_after!,
        ms_switching_condition,
        ms_bc, MS_p,
        (MS_ε, 10.0*MS_I0),
        1
    )

    nlsol = SwitchingControl.solve(
        ms_prob;
        mode = :free_boundary,
        guess = [0.25],
        make_u0 = ms_make_u0,
        residual! = ms_residual!
    )

    A_star = nlsol.u[1]
    u0_star = ms_make_u0([A_star])

    final_prob = SwitchingProblem(
        ms_dynamics_before!, ms_dynamics_after!, ms_switching_condition,
        BoundaryConditions(u0_star, nothing), MS_p, (MS_ε, 10.0*MS_I0), 1
    )
    final_sol = SwitchingControl.solve(final_prob)
    Vstar_numeric = final_sol.sol.tau

    rel_error = abs(Vstar_numeric - MS_Vstar_closed_form) / MS_Vstar_closed_form

    @test rel_error < 1e-4
end


# ============================================================
# DIXIT (1989) ENTRY/EXIT — MULTI-LEG FREE-BOUNDARY PROBLEM
#
# Two independently-integrated legs (idle region V0, active region
# V1) glued by value-matching + smooth-pasting at two endogenous
# thresholds (P_H, P_L). This is the only existing test that
# exercises `Leg`/`solve_multi_leg`, so it's kept in the main suite
# rather than examples/ despite being the most involved test here —
# regressions in the multi-leg machinery would otherwise go
# undetected. Tolerances are deliberately looser than the
# interactive/example version (which used 1e-12 ODE tolerances to
# show near-machine-precision agreement) purely to keep CI fast;
# 1e-8 still comfortably distinguishes a working solve from a
# broken one.
#
# Helper functions are prefixed `dx_` to avoid collisions with
# other testsets in this file.
# ============================================================

const DX_r = 0.05
const DX_δ = 0.02
const DX_σ = 0.20
const DX_w = 1.0
const DX_I = 0.7
const DX_E = 0.3

const DX_p = (r = DX_r, δ = DX_δ, σ = DX_σ, w = DX_w, I = DX_I, E = DX_E)

const DX_A2 = 0.5*DX_σ^2
const DX_A1 = (DX_r - DX_δ) - 0.5*DX_σ^2
const DX_A0 = -DX_r
const DX_disc = sqrt(DX_A1^2 - 4*DX_A2*DX_A0)
const DX_β1 = (-DX_A1 + DX_disc) / (2*DX_A2)
const DX_β2 = (-DX_A1 - DX_disc) / (2*DX_A2)

const DX_P_H_closed = 1.2701531078
const DX_P_L_closed = 0.7987448046

const DX_ε       = 0.10
const DX_P_upper = 5.0
const DX_P_far   = 50.0
const DX_P_lower = 0.01

function dx_dynamics_v0!(du, u, p, P)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*P*Fp) / (p.σ^2 * P^2)
    return nothing
end

function dx_dynamics_v1!(du, u, p, P)
    F, Fp = u[1], u[2]
    du[1] = Fp
    du[2] = (2*p.r*F - 2*(p.r - p.δ)*P*Fp - 2*(P - p.w)) / (p.σ^2 * P^2)
    return nothing
end

function dx_residual!(F, sols, θ, p)
    a, b, P_H, P_L = θ

    sol_v0, sol_v1 = sols[1], sols[2]

    V0_at_PH = sol_v0.sol.solution_before.u[end]
    V1_at_PL = sol_v1.sol.solution_before.u[end]

    V0_at_PL = sol_v0.sol.solution_before(P_L)
    V1_at_PH = sol_v1.sol.solution_before(P_H)

    F[1] = V1_at_PH[1] - p.I - V0_at_PH[1]
    F[2] = V1_at_PH[2] - V0_at_PH[2]
    F[3] = V0_at_PL[1] - (V1_at_PL[1] - p.E)
    F[4] = V0_at_PL[2] - V1_at_PL[2]

    return nothing
end

@testset "Dixit multi-leg free boundary" begin

    dx_v0_leg = Leg(
        dynamics! = dx_dynamics_v0!,
        tspan = (DX_ε, DX_P_upper),
        parameters = DX_p,
        make_u0 = θ -> begin
            a = θ[1]
            [a*DX_ε^DX_β1, a*DX_β1*DX_ε^(DX_β1 - 1)]
        end,
        switching_condition = θ -> begin
            P_H_guess = θ[3]
            (u, p, P) -> P - P_H_guess
        end
    )

    dx_v1_leg = Leg(
        dynamics! = dx_dynamics_v1!,
        tspan = (DX_P_far, DX_P_lower),
        parameters = DX_p,
        make_u0 = θ -> begin
            b = θ[2]
            [b*DX_P_far^DX_β2 + DX_P_far/DX_p.δ - DX_p.w/DX_p.r,
             b*DX_β2*DX_P_far^(DX_β2 - 1) + 1/DX_p.δ]
        end,
        switching_condition = θ -> begin
            P_L_guess = θ[4]
            (u, p, P) -> P - P_L_guess
        end
    )

    guess = [32.0, 2.8, 1.2, 0.85]

    nlsol = SwitchingControl.solve(
        [dx_v0_leg, dx_v1_leg], guess;
        residual! = dx_residual!,
        autodiff = :central,
        ivp_kwargs = (abstol = 1e-10, reltol = 1e-10),
        nlsolve_kwargs = (abstol = 1e-9, maxiters = 200)
    )

    a_star, b_star, P_H_star, P_L_star = nlsol.u

    # NonlinearSolve may return Stalled after reaching the numerical
    # solution because the finite-difference Jacobian is limited by the
    # inner ODE tolerance. Validate the actual free boundaries instead.
    @test all(isfinite, nlsol.u)

    @test abs(P_H_star - DX_P_H_closed) / DX_P_H_closed < 1e-8
    @test abs(P_L_star - DX_P_L_closed) / DX_P_L_closed < 1e-8
end

@testset "Guo & Zhang: σ1 sensitivity (Table 1)" begin
    σ1_grid  = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    x1_exact = [0.646, 0.531, 0.441, 0.369, 0.312, 0.266]
    x2_exact = [0.764, 0.683, 0.614, 0.554, 0.505, 0.462]

    for (σ1, x1e, x2e) in zip(σ1_grid, x1_exact, x2_exact)
        prob = MarkovSwitchingProblem(
            [GBMRegime(3.0, σ1), GBMRegime(3.0, 5.0)],
            [-100.0 100.0; 100.0 -100.0],
            [PutPayoff(5.0), PutPayoff(5.0)],
            3.0, (1e-3, 15.0), nothing, 2
        )
        sol = SwitchingControl.solve(prob; guess=[x1e*0.95, x2e*0.95])
        @test sol.converged
        @test isapprox(sol.boundaries[1], x1e; atol=0.015)
        @test isapprox(sol.boundaries[2], x2e; atol=0.015)
    end
end

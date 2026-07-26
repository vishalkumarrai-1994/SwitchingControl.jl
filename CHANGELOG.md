# Changelog

## v0.2.0
- Add `MarkovSwitchingProblem`: coupled N-regime free-boundary solver
  under a Markov generator, sharing the `solve()` interface via
  multiple dispatch on `AbstractSwitchingProblem`.
- Exact closed-form path (`GBMRegime` + `PutPayoff`, S=2) validated
  against Guo & Zhang (2004), base case + full σ1/λ1 sensitivity
  tables (24 parameter points).
- Generic N-regime / non-GBM path stubbed for upcoming numeric solver.

## v0.1.1 (untagged)
- Add `solve_free_boundary`: generalized shooting for endogenous
  free-boundary problems (smooth-pasting / value-matching).
- Add `solve_multi_leg`/`Leg` for multi-leg problems (Dixit entry/exit).
- Add McDonald-Siegel and Dixit worked examples, validated against
  closed-form solutions.

## v0.1.0
- Initial public release: ODE-based endogenous regime switching,
  ContinuousCallback event handling, IVP solver support.
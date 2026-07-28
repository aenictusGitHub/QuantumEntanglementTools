# Optimization and solvers

Optimization-backed APIs are not implemented at this scaffold milestone. ADR
0005 remains proposed until representative formulations are prototyped.

## Contract

Every optimization-based public result must expose:

- model/formulation identifier and mathematical direction;
- modeling layer, solver, and versions;
- termination plus primal and dual status;
- primal and dual objectives, gaps, and relevant residuals;
- tolerances, iteration count, time/budget, and stopping reason;
- certificate or bound interpretation;
- warnings and optional raw result.

`UNKNOWN`, time-limited, infeasible-or-unbounded, numerical failure, or
inaccurate statuses cannot be converted into a mathematical yes/no answer.

## Formulation gate

Before implementing an upstream CVX routine:

1. write its primal and, where available, dual formulation;
2. state cone and complex-variable requirements;
3. identify all convention conversions;
4. choose a package-owned high-level result;
5. test small primal/dual or independent formulations;
6. cover infeasible, unbounded, ill-conditioned, inaccurate, and time-limited
   outcomes.

Installed command-line GLPK is environment evidence only. It is not an SDP
solver and no solver backend is validated merely because an executable exists.

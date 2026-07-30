# Positive-semidefinite model constraints

QETLAB's `IsPSD` has two different behaviors: it tests ordinary numeric
matrices, but it adds a constraint when its input is a CVX expression.
QuantumEntanglementTools keeps those operations explicit:

- [`is_positive_semidefinite`](@ref) is the numeric tri-state predicate;
- [`positive_semidefinite_constraint`](@ref) creates a solver-neutral affine
  PSD constraint.

This separation prevents model construction from looking like a Boolean
mathematical conclusion.

## Constructing a constraint

The builder accepts a [`HermitianAffineMatrix`](@ref), copies its constant and
coefficient matrices, and returns a block that can be placed in
[`SemidefiniteProgram`](@ref). It does not create a JuMP model or select a
solver.

```julia
using QuantumEntanglementTools

rho = hermitian_variable(:rho, 2)
rho_is_psd = positive_semidefinite_constraint(rho)

trace_rho = trace_affine(rho)
trace_one = AffineEquality(
    AffineScalar(trace_rho.constant - 1, trace_rho.coefficients),
    :trace_one,
)
objective = AffineScalar(0.0, zeros(4))

problem = SemidefiniteProgram(
    :density_matrix_feasibility,
    :minimize,
    4,
    objective;
    equalities=[trace_one],
    psd_constraints=[rho_is_psd],
    primal_views=[rho_is_psd],
)
```

The optional JuMP extension translates that package-owned problem only after
the caller supplies an explicit backend:

```julia
using JuMP
using Hypatia

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)
result = solve_optimization(problem, backend)
```

When composing directly in an existing JuMP model, the extension also
provides the explicit materialization form
`positive_semidefinite_constraint(model, rho_is_psd, coordinates;
allow_densify=true)`. It returns a JuMP constraint reference and never calls
`optimize!`.

`allow_densify=true` is required because the extension represents a complex
Hermitian cone with a dense real symmetric block. Model limits are checked
before that allocation.

## Interpretation

`positive_semidefinite_constraint(rho)` asserts no fact about a particular
numeric matrix. It only returns model data. Likewise, a solver's
`OptimizationOptimal` status is not automatically a theorem-level
certificate: inspect its objective, residuals, tolerances, backend metadata,
and reconstructed primal and dual evidence.

For a numeric array `A`, use:

```julia
result = is_positive_semidefinite(A)
```

That predicate never normalizes, Hermitian-symmetrizes, projects, clips, or
repairs `A`. Inputs outside its validation or resource limits remain explicit
errors or `unknown`, not negative mathematical results.

## Upstream correspondence

This is the Julia-native replacement for the CVX-expression branch of QETLAB
`IsPSD.m` at the pinned upstream revision. The numeric branch remains mapped to
`is_positive_semidefinite`. The Julia API intentionally does not reproduce
QETLAB's ambient-model mutation: constraint construction is explicit and the
returned matrices are owned by the package model.

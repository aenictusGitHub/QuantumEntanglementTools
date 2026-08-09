# Executable tutorials

These tutorials are ordinary Julia scripts from the repository's `tutorials/`
directory. Each script performs live calculations, checks its mathematical
invariants, prints values derived from that run, and returns a small summary for
the automated tutorial gate.

## Choose a tutorial

The suggested order moves from array and subsystem conventions to
certificate-aware entanglement examples. Runtime estimates describe the
calculation after Julia and the package have loaded; first-time compilation
depends on the machine.

| Order | Tutorial | Level | Typical runtime | Optional dependencies |
|---:|---|---|---|---|
| 1 | Subsystem reductions and ordering | Beginner | Seconds | None |
| 2 | Seeded Schmidt decomposition | Beginner | Seconds | None |
| 3 | Local channel noise | Intermediate | Seconds | None |
| 4 | Separability certificates | Beginner | Seconds | None |
| 5 | Certificates and `unknown` | Intermediate | Seconds | None |
| 6 | Tiles PPT bound entanglement | Advanced | Seconds | None |
| 7 | Symmetric SAPPT witnesses | Advanced | Seconds | None |

Run any tutorial from the repository root:

```sh
julia --startup-file=no --project=. tutorials/subsystem_reductions.jl
julia --startup-file=no --project=. tutorials/qetlab_intro_schmidt.jl
julia --startup-file=no --project=. tutorials/local_channel_noise.jl
julia --startup-file=no --project=. tutorials/separability_examples.jl
julia --startup-file=no --project=. tutorials/entanglement_certificates.jl
julia --startup-file=no --project=. tutorials/qetlab_intro_tiles.jl
julia --startup-file=no --project=. tutorials/symmetric_sappt_witnesses.jl
```

Run the complete tutorial gate with:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

The package test entry point includes that same runner, so `Pkg.test()` also
executes every tutorial. None of these workflows needs an optional backend.
The short examples below expose the essential API calls; the `include(...)`
blocks execute the complete assertion-backed scripts.

## 1. Subsystem reductions and ordering

The first workflow constructs $|0\rangle_A \otimes |\Phi^+\rangle_{BC}$.
It traces out selected systems, moves the product qubit from the first position
to the last, and checks the Bell pair's partial-transpose spectrum. This makes
the package's left-to-right subsystem order observable rather than implicit.

A minimal reduction is directly copyable:

```@example tutorial-subsystems-core
using QuantumEntanglementTools

psi = tensor_product([1.0, 0.0], bell_state())
rho_BC = partial_trace(psi, (2, 2, 2); trace_out=(1,))
size(rho_BC), purity(rho_BC)
```

Run the complete tutorial:

```@example tutorial-subsystems
using QuantumEntanglementTools

path = joinpath(pkgdir(QuantumEntanglementTools), "tutorials", "subsystem_reductions.jl")
include(path);
TutorialSubsystemReductions.run()
```

The vector overload of `partial_trace` forms the mathematical pure-state
reduction without requiring the tutorial to materialize the full three-qubit
density matrix. See [Mathematical conventions](conventions.md) for subsystem
and basis ordering.

## 2. A seeded 3 × 3 Schmidt decomposition

This independently written workflow follows the mathematical arc of the
[QETLAB homepage](https://qetlab.com/) introduction. A fixed `Xoshiro` seed
produces a reproducible random pure state on
$\mathbb{C}^3 \otimes \mathbb{C}^3$. The tutorial computes its Schmidt
decomposition, reconstructs the state both term by term and with
`tensor_sum(...; weights=...)`, and checks normalization, orthonormal Schmidt
vectors, rank, and reconstruction residuals.

The essential calls are:

```@example tutorial-qetlab-schmidt-core
using QuantumEntanglementTools
using Random: Xoshiro

psi = random_state_vector(Xoshiro(2026), (3, 3))
decomposition = schmidt_decomposition(psi, (3, 3))
decomposition.coefficients
```

Run the complete tutorial:

```@example tutorial-qetlab-schmidt
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools), "tutorials", "qetlab_intro_schmidt.jl"
)
include(path);
TutorialQETLABIntroSchmidt.run()
```

The explicit random-number generator makes the example repeatable without
changing Julia's global random stream. The workflow was written independently
from the documented contracts of QETLAB's `RandomStateVector.m`,
`SchmidtDecomposition.m`, `Tensor.m`, and `TensorSum.m` at pinned revision
`d8589610f00cff106537268dee2e2a1153f3a601`; attribution and BSD-2-Clause
license details are recorded in `PROVENANCE.toml` and
`licenses/QETLAB-LICENSE.txt`.

## 3. A local channel acting on an entangled state

This workflow creates a completely depolarizing qubit channel, verifies that it
is completely positive, trace preserving, and unital, and compares its Choi and
superoperator representations. Applying it to one half of a Bell state produces
the maximally mixed two-qubit state.

```@example tutorial-channel
using QuantumEntanglementTools

path = joinpath(pkgdir(QuantumEntanglementTools), "tutorials", "local_channel_noise.jl")
include(path);
TutorialLocalChannelNoise.run()
```

Both entanglement conclusions in this workflow name their certificates. The
initial negative-partial-transpose witness certifies entanglement; after local
depolarization, the PPT theorem in the exact $2\times2$ domain certifies
separability.

## 4. Separability certificates by example

This workflow starts with a visibly factorized pure state, builds full-rank
$2\times2$ and $3\times3$ density matrices from product-basis projectors,
and compares the native criterion pipeline with the sufficient separable-ball
test. It also shows why `:outside_ball` is not an entanglement verdict.

The shortest certificate-bearing example is:

```@example tutorial-separability-core
using QuantumEntanglementTools

psi = tensor_product(ComplexF64[1, 0], ComplexF64[0, 1])
analyze_entanglement(psi, (2, 2))
```

Run the complete tutorial:

```@example tutorial-separability
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools), "tutorials", "separability_examples.jl"
)
include(path);
TutorialSeparabilityExamples.run()
```

The examples deliberately include one state for which
`analyze_entanglement` returns `:unknown` while `in_separable_ball` returns
`:separable_certified`. See [Separability by example](separability_examples.md)
for the construction line by line and a guide to the two result vocabularies.

## 5. Certificates, necessary tests, and `unknown`

The certificate workflow contrasts four outcomes:

- Schmidt rank certifies that a pure Bell state is entangled.
- An explicit product decomposition certifies a pure product state as
  separable.
- The $3\times3$ Horodecki example passes the PPT attempt without obtaining a
  separability certificate, then violates the realignment criterion.
- The maximally mixed $3\times3$ state passes all requested necessary tests,
  but the native pipeline conservatively reports `unknown` because those passes
  do not constitute a separability certificate in that dimension.

```@example tutorial-certificates
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools), "tutorials", "entanglement_certificates.jl"
)
include(path);
TutorialEntanglementCertificates.run()
```

Inspect `report.attempts` whenever the route to a conclusion matters. An
`unknown` result is a deliberate statement about available certification, not
a synonym for separable or entangled. See [Entanglement
backends](entanglement_backends.md) for the complete result semantics.

## 6. Tiles: an exact PPT bound-entanglement proof

The Tiles workflow starts from five product vectors with integer local factors.
It verifies their unextendibility exactly, constructs the rational projector
onto their span, and normalizes the complementary rank-four projector to obtain
a density matrix. Exact arithmetic then proves that the density matrix is fixed
by partial transpose. Its range is the complementary subspace, while the exact
UPB certificate says that this subspace contains no product vector; the range
criterion therefore proves entanglement. Together, these two facts prove that
the state is PPT bound entangled.

```@example tutorial-qetlab-tiles
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools), "tutorials", "qetlab_intro_tiles.jl"
)
include(path);
TutorialQETLABIntroTiles.run()
```

The floating-point checks deliberately preserve the distinction between an
exact proof and a numerical test. `ppt_criterion` returns `CriterionUnknown`
for the rank-four state because it lies on the zero-eigenvalue boundary; this
does not weaken the exact partial-transpose identity. For an independent,
platform-stable numerical check, the tutorial mixes in exactly `1//1024` of the
maximally mixed state before converting to floating point. This full-rank
neighbor remains PPT by the same exact partial-transpose identity. The combined
`is_separable(...; strategies=(:ppt, :realignment))` call records PPT as a
necessary but insufficient test and then certifies the neighbor's entanglement
through a realignment cross-norm violation. Because the maximally mixed state
is separable and separable states form a convex set, entanglement of this
depolarized neighbor also independently implies entanglement of the original
boundary state.

This workflow was also written independently from the Tiles example on the
[QETLAB homepage](https://qetlab.com/) and the documented contracts of
`UPB.m`, `IsUPB.m`, `IsPPT.m`, `IsSeparable.m`, `Tensor.m`, and
`PartialTranspose.m` at the pinned revision above. The same repository
provenance and license records apply.

## 7. Symmetric SAPPT states and witnesses

This workflow reproduces the five-qubit family studied by Louvet *et al.* It
constructs a 19-term separable decomposition, compares it with a same-spectrum
GHZ representative, reconstructs the published symmetric witness, and builds a
decomposable NPT witness below the SAPPT threshold. The code also treats the
GHZ phase convention explicitly.

```@example tutorial-symmetric-sappt
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools), "tutorials", "symmetric_sappt_witnesses.jl"
)
include(path);
TutorialSymmetricSAPPTWitnesses.run()
```

The rounded witness coefficients reproduce the paper's reported values; the
optimization that originally produced them is not rerun. See
[Symmetric SAPPT states and witnesses](paper_symmetric_separability.md) for the
state family, finite separable decomposition, block-positivity calculation,
phase convention, and certification boundaries.

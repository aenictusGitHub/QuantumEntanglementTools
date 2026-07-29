# Executable tutorials

These tutorials are ordinary Julia scripts from the repository's `tutorials/`
directory. Each script performs live calculations, checks its mathematical
invariants, prints values derived from that run, and returns a small summary for
the automated tutorial gate.

Run any tutorial from the repository root:

```sh
julia --startup-file=no --project=. tutorials/subsystem_reductions.jl
julia --startup-file=no --project=. tutorials/local_channel_noise.jl
julia --startup-file=no --project=. tutorials/entanglement_certificates.jl
julia --startup-file=no --project=. tutorials/separability_examples.jl
julia --startup-file=no --project=. tutorials/symmetric_sappt_witnesses.jl
```

Run the complete tutorial gate with:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

The package test entry point includes that same runner, so `Pkg.test()` also
executes every tutorial. None of these workflows needs an optional backend.

## Subsystem reductions and ordering

The first workflow constructs $|0\rangle_A \otimes |\Phi^+\rangle_{BC}$.
It traces out selected systems, moves the product qubit from the first position
to the last, and checks the Bell pair's partial-transpose spectrum. This makes
the package's left-to-right subsystem order observable rather than implicit.

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

## A local channel acting on an entangled state

The second workflow creates a completely depolarizing qubit channel, verifies
that it is completely positive, trace preserving, and unital, and compares its
Choi and superoperator representations. Applying it to one half of a Bell state
produces the maximally mixed two-qubit state.

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

## Separability certificates by example

This workflow starts with a visibly factorized pure state, builds full-rank
$2\times2$ and $3\times3$ density matrices from product-basis projectors,
and compares the native criterion pipeline with the sufficient separable-ball
test. It also shows why `:outside_ball` is not an entanglement verdict.

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

## Symmetric SAPPT states and witnesses

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

## Certificates, necessary tests, and `unknown`

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

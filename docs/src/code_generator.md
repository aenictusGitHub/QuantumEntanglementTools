# Entanglement example code generator

Build a complete Julia example from a curated state family and a set of
certificate-aware analyses. The generator accepts only bounded numeric values
and fixed choices; it does not accept or evaluate arbitrary Julia code.

> **Viewing this source on GitHub?** GitHub renders the Markdown but does not
> execute the Documenter assets that power the interactive form. The generator
> runs in the
> [live rendered documentation](https://aenictusgithub.github.io/QuantumEntanglementTools/code_generator/).
> If Pages has not yet been enabled for a new checkout, open the
> [Documentation workflow](https://github.com/aenictusGitHub/QuantumEntanglementTools/actions/workflows/docs.yml),
> download and unzip its
> `documentation-<run-id>-<run-attempt>` artifact, then open
> `code_generator/index.html`. Artifacts are retained for seven days.
>
> To build the same documentation from a local checkout, run:
>
> ```sh
> julia --startup-file=no scripts/build_docs.jl
> ```
>
> The helper uses a temporary environment for the running Julia version and
> exits with an error if resolution or the strict build fails. Then open
> `docs/build/code_generator.html` in a browser. For examples that run without a
> browser, use the [executable tutorials](tutorials.md).

## Interactive generator

If the form is not visible, continue to the [static fallback](#Static-fallback)
below or use an [executable tutorial](tutorials.md); the Markdown source remains
useful even when browser assets are unavailable.

## Static fallback

If the interactive assets are unavailable, this complete example exercises the
same certificate-aware Bell-state path and can be pasted directly into Julia:

```julia
using QuantumEntanglementTools

psi = bell_state()
rho = psi * psi'
dims = (2, 2)

report = analyze_entanglement(rho, dims; atol=1e-12, rtol=0)
ppt = ppt_criterion(rho, dims; atol=1e-12, rtol=0)

@assert conclusion(report) === :entangled
@assert is_certified(report)
@assert ppt.status === CriterionEntanglementDetected

println(explain(report))
```

The executable version does not upload data or evaluate dynamically supplied
code. Additional complete scripts are listed under
[Executable tutorials](tutorials.md).

## Generator capabilities

The curated state menu includes pure product, Bell, GHZ, and Dicke vectors;
explicit diagonal, isotropic, Werner, Horodecki, and five-qubit symmetric mixed
states; the exact 3 × 3 Tiles-UPB construction; and seeded random density
matrices. Random examples use a fixed local seed, bounded dimensions and rank,
a Hilbert–Schmidt or Bures choice, and an optional real-valued ensemble. They
never consume Julia's caller-global random stream.

The analysis menu separates four kinds of work:

- non-mutating density validation, marginals, and pure-state Schmidt data;
- certificate-aware pipeline and structured separability reports;
- individual PPT, realignment, reduction, separable-ball, and guarded-witness
  routes;
- measures and side-effect-free optional-backend readiness.

Schmidt output is available only for curated pure-state families. Separability
profiles expand to fixed deterministic core strategies; selecting a profile
cannot inject a function name or silently load an optional solver. Output
options can use the package's concise result helpers and add bounded structural
assertions to the generated script.

The resource card is a live planning estimate of the chosen bipartition,
Hilbert dimension, dense entry count, and approximate matrix storage.
Generation remains the authoritative validation step.

## Saving and loading a configuration

**Copy configuration** and **Download configuration** export canonical,
versioned JSON through the same normalization used for code generation.
**Load a JSON configuration** accepts at most 16 KiB and passes its contents to
the generator's strict parser before changing any control. A loaded
configuration is not executed and does not regenerate Julia automatically:
review it, then press **Generate Julia code**.

Configuration schema 2 records the family, bounded parameters, analyses,
strategy profile, tolerances, and output choices. Schema-1 configurations are
accepted with compatible defaults for the new options. Unknown future schema
versions, unknown keys or choices, malformed JSON, non-finite numbers, and
out-of-range values are rejected. The file name is fixed to
`qet_generator_config.json`; configurations are not placed in a URL, sent to a
server, evaluated as code, or saved in browser storage.

## Interpreting the output

The generated program preserves the package's three-valued and certificate-aware
semantics:

| Output | What it establishes |
|---|---|
| Pipeline `certified == true` | The named certificate supports the reported `:separable` or `:entangled` conclusion. |
| Separability report with `is_certified(result) == true` | The named strategy produced a certificate; inspect `conclusion(result)` and `explain(result)` for its scope. |
| `CriterionEntanglementDetected` | A robust PPT, realignment, or reduction violation certifies entanglement. |
| `CriterionSatisfied` | A necessary condition passed; this is not a general separability conclusion. |
| Ball status `:separable_certified` | The sufficient separable-ball inequality certifies separability. |
| `:unknown` or `:outside_ball` | This computation is inconclusive. |

Density validation reports defects but never repairs, normalizes, clips, or
symmetrizes the state. Backend readiness reports availability only; it does not
select a backend, solve an optimization problem, or establish separability.

The decomposable PPT witness is emitted only inside a
`CriterionEntanglementDetected` guard. The five-qubit option separately
constructs the rounded symmetric-subspace witness documented in
[Symmetric SAPPT states and witnesses](paper_symmetric_separability.md). A
negative expectation is conclusive; a nonnegative expectation is not.

## Running a downloaded example

The package is not registered. Run a downloaded file from an environment in
which this checkout has been developed:

```julia
using Pkg
Pkg.develop(path="/path/to/QuantumEntanglementTools")
```

Then use the command printed at the top of the generated file, for example:

```sh
julia --startup-file=no --project=. qet_bell_example.jl
```

Generator limits—local dimensions at most eight, general dense Hilbert
dimension at most 256, random-state Hilbert dimension at most 64, rank no
larger than the selected random state's dimension, structured separability
reports capped at Hilbert dimension 64, a bounded integer seed, the fixed
3 × 3 Tiles construction, and the fixed five-qubit paper family—bound the cost
of emitted spectral workflows. Configuration files are additionally capped at
16 KiB. These are safety limits of this documentation tool, not claims about
the package's mathematical API limits.

## License of generated artifacts

The emitted Julia program contains substantial generator-supplied template
text. It is covered by this repository's **BSD-3-Clause** license and carries
`SPDX-FileCopyrightText` and `SPDX-License-Identifier` headers plus a link to
the complete license. Retain those notices and satisfy the BSD terms when
redistributing a generated program or substantial portions of its template.

The downloaded JSON configuration is descriptive metadata containing
normalized user selections; it does not embed the Julia or JavaScript
template. The generator does not claim rights in parameters, comments, data,
or other material supplied by a user. Users remain responsible for having the
right to distribute material they add to either artifact.

For fully prewritten and continuously tested alternatives, see
[Executable tutorials](tutorials.md) and
[Separability by example](separability_examples.md).

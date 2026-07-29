# Entanglement example code generator

Build a complete Julia example from a curated state family and a set of
certificate-aware analyses. The generator accepts only bounded numeric values
and fixed choices; it does not accept or evaluate arbitrary Julia code.

> **Viewing this source on GitHub?** GitHub renders the Markdown but does not
> execute the Documenter assets that power the interactive form. The generator
> is available only in built documentation. For a recent pushed commit, open
> the [Documentation workflow](https://github.com/aenictusGitHub/QuantumEntanglementTools/actions/workflows/docs.yml),
> download and unzip its `documentation` artifact, then open
> `code_generator/index.html`. Artifacts are retained for seven days.
>
> To build the same documentation from a local checkout, run:
>
> ```sh
> julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
> julia --startup-file=no --project=docs docs/make.jl
> ```
>
> Then open `docs/build/code_generator.html` in a browser. For examples that run
> without a browser, use the [executable tutorials](tutorials.md).

## Interactive generator

## Interpreting the output

The generated program preserves the package's three-valued and certificate-aware
semantics:

| Output | What it establishes |
|---|---|
| Pipeline `certified == true` | The named certificate supports the reported `:separable` or `:entangled` conclusion. |
| `CriterionEntanglementDetected` | A robust PPT, realignment, or reduction violation certifies entanglement. |
| `CriterionSatisfied` | A necessary condition passed; this is not a general separability conclusion. |
| Ball status `:separable_certified` | The sufficient separable-ball inequality certifies separability. |
| `:unknown` or `:outside_ball` | This computation is inconclusive. |

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

Generator limits—local dimensions at most eight, GHZ/Dicke dense Hilbert
dimension at most 256, and the fixed five-qubit paper family—bound the cost of
the emitted spectral workflows. They are safety limits of this documentation
tool, not claims about the package's mathematical API limits.

The original generator templates are part of this BSD-3-Clause repository.
Retain the repository's license notice when redistributing substantial template
text. For fully prewritten and continuously tested alternatives, see
[Executable tutorials](tutorials.md) and
[Separability by example](separability_examples.md).

# Entanglement example code generator

Build a complete Julia example from a curated state family and a set of
certificate-aware analyses. The generator accepts only bounded numeric values
and fixed choices; it does not accept or evaluate arbitrary Julia code.

```@raw html
<div id="qet-code-generator" class="qet-generator" data-generator-schema="1">
  <div class="qet-generator__local-note" role="note">
    <strong>Local and deterministic.</strong>
    This page generates text in your browser. It does not run Julia, send the
    selected parameters to a server, or save them after the page is reloaded.
  </div>

  <noscript>
    <div class="qet-generator__noscript" role="alert">
      JavaScript is required for the interactive builder. The executable
      tutorials linked below remain available without it.
    </div>
  </noscript>

  <div class="qet-generator__layout">
    <form id="qet-generator-form" class="qet-generator__form" novalidate>
      <fieldset class="qet-generator__section">
        <legend>1. Start from a tested preset</legend>
        <div class="qet-generator__field">
          <label for="qet-generator-preset">Preset</label>
          <select id="qet-generator-preset"></select>
          <p id="qet-preset-description" class="qet-generator__hint"></p>
        </div>
      </fieldset>

      <fieldset class="qet-generator__section">
        <legend>2. Construct the state</legend>
        <div class="qet-generator__field">
          <label for="qet-state-family">State family</label>
          <select id="qet-state-family"></select>
          <p class="qet-generator__field-error" data-qet-error-for="qet-state-family" hidden></p>
        </div>

        <div class="qet-generator__family" data-qet-family="product_basis">
          <p class="qet-generator__hint">
            Construct <span aria-label="ket i sub A tensor ket j sub B">|i<sub>A</sub>⟩ ⊗ |j<sub>B</sub>⟩</span>.
            Basis labels below are zero-based, while the emitted Julia array
            positions are one-based.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-product-dim-a">Dimension of A</label>
              <input id="qet-product-dim-a" type="number" value="2" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-product-dim-a" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-product-index-a">Basis label i</label>
              <input id="qet-product-index-a" type="number" value="0" min="0" max="1" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-product-index-a" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-product-dim-b">Dimension of B</label>
              <input id="qet-product-dim-b" type="number" value="2" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-product-dim-b" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-product-index-b">Basis label j</label>
              <input id="qet-product-index-b" type="number" value="1" min="0" max="1" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-product-index-b" hidden></p>
            </div>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="bell" hidden>
          <div class="qet-generator__field">
            <label for="qet-bell-index">Bell-state index</label>
            <select id="qet-bell-index">
              <option value="0">0 — (|00⟩ + |11⟩) / √2</option>
              <option value="1">1 — (|00⟩ − |11⟩) / √2</option>
              <option value="2">2 — (|01⟩ + |10⟩) / √2</option>
              <option value="3">3 — (|01⟩ − |10⟩) / √2</option>
            </select>
            <p class="qet-generator__field-error" data-qet-error-for="qet-bell-index" hidden></p>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="diagonal_mixture" hidden>
          <p class="qet-generator__hint">
            Emit every product-basis vector, its projector, and a binary-exact
            convex weight. The 3 × 3 choice illustrates an inconclusive native
            pipeline together with a sufficient separable-ball certificate.
          </p>
          <div class="qet-generator__field">
            <label for="qet-mixture-dimension">Local dimension</label>
            <select id="qet-mixture-dimension">
              <option value="2">2 × 2</option>
              <option value="3" selected>3 × 3</option>
            </select>
            <p class="qet-generator__field-error" data-qet-error-for="qet-mixture-dimension" hidden></p>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="ghz" hidden>
          <p class="qet-generator__hint">
            The emitted bipartition groups the first <em>k</em> parties against
            the remaining parties without changing the package's basis order.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-ghz-dimension">Local dimension</label>
              <input id="qet-ghz-dimension" type="number" value="2" min="2" max="4" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-ghz-dimension" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-ghz-parties">Number of parties</label>
              <input id="qet-ghz-parties" type="number" value="3" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-ghz-parties" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-ghz-cut">Parties left of cut</label>
              <input id="qet-ghz-cut" type="number" value="1" min="1" max="2" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-ghz-cut" hidden></p>
            </div>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="dicke" hidden>
          <p class="qet-generator__hint">
            Construct a normalized qubit Dicke state and group a contiguous
            left block against the remaining qubits.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-dicke-parties">Number of qubits</label>
              <input id="qet-dicke-parties" type="number" value="4" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-dicke-parties" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-dicke-excitations">Excitations</label>
              <input id="qet-dicke-excitations" type="number" value="2" min="0" max="4" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-dicke-excitations" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-dicke-cut">Qubits left of cut</label>
              <input id="qet-dicke-cut" type="number" value="2" min="1" max="3" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-dicke-cut" hidden></p>
            </div>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="isotropic" hidden>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-isotropic-dimension">Local dimension</label>
              <input id="qet-isotropic-dimension" type="number" value="3" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-isotropic-dimension" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-isotropic-alpha">Mixing parameter α</label>
              <input id="qet-isotropic-alpha" type="number" value="0.5" min="-1" max="1" step="any">
              <p class="qet-generator__field-error" data-qet-error-for="qet-isotropic-alpha" hidden></p>
            </div>
          </div>
          <p class="qet-generator__hint">
            The exact positivity lower bound −1/(d²−1) is checked when code is generated.
          </p>
        </div>

        <div class="qet-generator__family" data-qet-family="werner" hidden>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-werner-dimension">Local dimension</label>
              <input id="qet-werner-dimension" type="number" value="3" min="2" max="8" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-werner-dimension" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-werner-alpha">Werner parameter α</label>
              <input id="qet-werner-alpha" type="number" value="0.8" min="-1" max="1" step="any">
              <p class="qet-generator__field-error" data-qet-error-for="qet-werner-alpha" hidden></p>
            </div>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="horodecki" hidden>
          <p class="qet-generator__hint">
            Interior 3 × 3 states are PPT entangled. Passing PPT is not the
            certificate; the default native pipeline detects the selected
            example through realignment.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-horodecki-a">Family parameter a</label>
              <input id="qet-horodecki-a" type="number" value="0.3" min="0" max="1" step="any">
              <p class="qet-generator__field-error" data-qet-error-for="qet-horodecki-a" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-horodecki-dims">Local dimensions</label>
              <select id="qet-horodecki-dims">
                <option value="3x3">3 × 3</option>
                <option value="2x4">2 × 4</option>
              </select>
              <p class="qet-generator__field-error" data-qet-error-for="qet-horodecki-dims" hidden></p>
            </div>
          </div>
        </div>

        <div class="qet-generator__family" data-qet-family="symmetric_sappt_ghz5" hidden>
          <p class="qet-generator__hint">
            Construct the five-qubit GHZ representative in the Dicke basis,
            embed it into the full Hilbert space, and build the paper's rounded
            phase-matched symmetric witness. The default p = 121/125 is a
            strict SAPPT point with negative witness expectation.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-symmetric-p">Mixing weight p</label>
              <input id="qet-symmetric-p" type="number" value="0.968" min="0" max="1" step="any">
              <p class="qet-generator__field-error" data-qet-error-for="qet-symmetric-p" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-symmetric-phase">GHZ phase (radians)</label>
              <input id="qet-symmetric-phase" type="number" value="0" step="any">
              <p class="qet-generator__field-error" data-qet-error-for="qet-symmetric-phase" hidden></p>
            </div>
          </div>
        </div>
      </fieldset>

      <fieldset class="qet-generator__section">
        <legend>3. Select analyses</legend>
        <p class="qet-generator__hint">
          Every spectral route uses a bounded dense state. No method is silently
          substituted for another.
        </p>
        <div class="qet-generator__checks">
          <label class="qet-generator__check">
            <input id="qet-analysis-pipeline" type="checkbox" checked>
            <span>Certificate-aware pipeline</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-pipeline" hidden></p>

          <label class="qet-generator__check">
            <input id="qet-analysis-ppt" type="checkbox">
            <span>PPT criterion</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-ppt" hidden></p>

          <label class="qet-generator__check qet-generator__check--nested">
            <input id="qet-analysis-ppt-witness" type="checkbox" disabled>
            <span>Construct decomposable PPT witness when negativity is robust</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-ppt-witness" hidden></p>

          <label class="qet-generator__check">
            <input id="qet-analysis-realignment" type="checkbox">
            <span>Realignment / CCNR criterion</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-realignment" hidden></p>

          <label class="qet-generator__check">
            <input id="qet-analysis-reduction" type="checkbox">
            <span>Reduction criterion on both sides</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-reduction" hidden></p>

          <label class="qet-generator__check">
            <input id="qet-analysis-separable-ball" type="checkbox">
            <span>Gurvits–Barnum separable-ball certificate</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-separable-ball" hidden></p>

          <label class="qet-generator__check">
            <input id="qet-analysis-measures" type="checkbox" checked>
            <span>Purity, entropy, negativity, and logarithmic negativity</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-measures" hidden></p>
        </div>

        <div class="qet-generator__field-grid qet-generator__tolerances">
          <div class="qet-generator__field">
            <label for="qet-analysis-atol">Absolute tolerance</label>
            <input id="qet-analysis-atol" type="text" inputmode="decimal" value="1e-12">
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-atol" hidden></p>
          </div>
          <div class="qet-generator__field">
            <label for="qet-analysis-rtol">Relative tolerance</label>
            <input id="qet-analysis-rtol" type="text" inputmode="decimal" value="0">
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-rtol" hidden></p>
          </div>
        </div>
        <p class="qet-generator__hint">
          Inputs are never normalized, clipped, symmetrized, or otherwise repaired.
        </p>
      </fieldset>

      <div id="qet-generator-errors" class="qet-generator__errors" role="alert" hidden></div>

      <div class="qet-generator__form-actions">
        <button class="qet-generator__button qet-generator__button--primary" type="submit">
          Generate Julia code
        </button>
        <button class="qet-generator__button" type="reset">Reset</button>
      </div>
    </form>

    <section class="qet-generator__output" aria-labelledby="qet-generator-output-title">
      <div class="qet-generator__output-heading">
        <div>
          <p class="qet-generator__eyebrow">Review and run</p>
          <h2 id="qet-generator-output-title">Generated Julia example</h2>
        </div>
        <div class="qet-generator__output-actions">
          <button id="qet-generator-copy" class="qet-generator__button" type="button" disabled>
            Copy code
          </button>
          <button id="qet-generator-download" class="qet-generator__button" type="button" disabled>
            Download .jl
          </button>
        </div>
      </div>

      <div id="qet-generator-summary" class="qet-generator__summary" aria-label="Generated example summary"></div>
      <div id="qet-generator-notices" class="qet-generator__notices" role="note" hidden></div>
      <p id="qet-generator-status" class="qet-generator__status" aria-live="polite">
        Loading the example builder…
      </p>
      <pre class="qet-generator__code" tabindex="0" aria-label="Generated Julia code"><code id="qet-generator-code"># Loading generator assets…
</code></pre>
    </section>
  </div>
</div>
```

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

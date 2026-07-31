(function () {
    "use strict";

    var generator = null;
    var root = null;
    var form = null;
    var lastResult = null;
    var applyingConfig = false;
    var CONFIG_MAX_BYTES = 16 * 1024;
    var CONFIG_FILENAME = "qet_generator_config.json";

    // Keep the full widget markup out of the Markdown source: GitHub does not
    // interpret Documenter's `@raw html` blocks, while built docs load this asset.
    var GENERATOR_MARKUP = `
<div id="qet-code-generator" class="qet-generator" data-generator-schema="2">
  <div class="qet-generator__local-note" role="note">
    <strong>Local and deterministic.</strong>
    This page generates text in your browser. It does not run Julia or send the
    selected parameters to a server. Nothing is persisted unless you explicitly
    copy or download generated code or a configuration.
  </div>


  <div class="qet-generator__layout">
    <form id="qet-generator-form" class="qet-generator__form" novalidate>
      <fieldset class="qet-generator__section">
        <legend>1. Start from a tested preset</legend>
        <div class="qet-generator__field">
          <label for="qet-generator-preset">Preset</label>
          <select id="qet-generator-preset" aria-describedby="qet-preset-description"></select>
          <p id="qet-preset-description" class="qet-generator__hint"></p>
        </div>
      </fieldset>

      <fieldset class="qet-generator__section">
        <legend>2. Construct the state</legend>
        <div class="qet-generator__field">
          <label for="qet-state-family">State family</label>
          <select id="qet-state-family" aria-describedby="qet-state-family-hint"></select>
          <p id="qet-state-family-hint" class="qet-generator__hint">
            Choose a curated construction. Parameters are checked against the
            generator's numeric and memory limits before code is emitted.
          </p>
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
            Interior states in the default 3 × 3 branch are PPT entangled.
            Passing PPT is not the certificate; the default native pipeline
            detects the selected example through realignment.
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

        <div class="qet-generator__family" data-qet-family="tiles_upb" hidden>
          <p class="qet-generator__hint">
            Build the exact 3 × 3 Tiles unextendible-product-basis construction
            and its normalized complementary state. The generated example keeps
            its PPT evidence separate from the range-criterion entanglement
            certificate.
          </p>
          <p class="qet-generator__fixed-choice">
            Fixed bipartition: <strong>3 × 3</strong>; dense Hilbert dimension:
            <strong>9</strong>.
          </p>
        </div>

        <div class="qet-generator__family" data-qet-family="random_density" hidden>
          <p id="qet-random-density-hint" class="qet-generator__hint">
            Construct a reproducible random density matrix from a fixed seed.
            Hilbert–Schmidt and Bures sampling are fixed, bounded choices; no
            user-supplied code or global random stream is used.
          </p>
          <div class="qet-generator__field-grid">
            <div class="qet-generator__field">
              <label for="qet-random-dim-a">Dimension of A</label>
              <input id="qet-random-dim-a" type="number" value="2" min="2" max="8" step="1" aria-describedby="qet-random-density-hint">
              <p class="qet-generator__field-error" data-qet-error-for="qet-random-dim-a" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-random-dim-b">Dimension of B</label>
              <input id="qet-random-dim-b" type="number" value="2" min="2" max="8" step="1" aria-describedby="qet-random-density-hint">
              <p class="qet-generator__field-error" data-qet-error-for="qet-random-dim-b" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-random-rank">Requested rank</label>
              <input id="qet-random-rank" type="number" value="4" min="1" max="4" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-random-rank" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-random-seed">Seed</label>
              <input id="qet-random-seed" type="number" value="2025" min="0" max="4294967295" step="1">
              <p class="qet-generator__field-error" data-qet-error-for="qet-random-seed" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-random-distribution">Distribution</label>
              <select id="qet-random-distribution">
                <option value="hilbert_schmidt">Hilbert–Schmidt</option>
                <option value="bures">Bures</option>
              </select>
              <p class="qet-generator__field-error" data-qet-error-for="qet-random-distribution" hidden></p>
            </div>
          </div>
          <label class="qet-generator__check qet-generator__check--standalone">
            <input id="qet-random-real" type="checkbox">
            <span>Use a real-valued ensemble</span>
          </label>
          <p class="qet-generator__field-error" data-qet-error-for="qet-random-real" hidden></p>
        </div>

        <aside id="qet-generator-resource-plan" class="qet-generator__resource-plan" aria-label="Resource planning estimate">
          Estimating the selected state's dense resource footprint…
        </aside>
      </fieldset>

      <fieldset class="qet-generator__section">
        <legend>3. Select analyses</legend>
        <p class="qet-generator__hint">
          Every spectral route uses a bounded dense state. No method is silently
          substituted for another.
        </p>
        <div class="qet-generator__analysis-group" role="group" aria-labelledby="qet-analysis-input-title">
          <h3 id="qet-analysis-input-title">Understand the input</h3>
          <div class="qet-generator__checks">
            <label class="qet-generator__check">
              <input id="qet-analysis-validation" type="checkbox">
              <span>Validate the density matrix without repairing it</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-validation" hidden></p>

            <label class="qet-generator__check">
              <input id="qet-analysis-marginals" type="checkbox">
              <span>Compute both reduced density matrices</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-marginals" hidden></p>

            <label class="qet-generator__check">
              <input id="qet-analysis-schmidt" type="checkbox" aria-describedby="qet-analysis-schmidt-hint">
              <span>Show the Schmidt decomposition (pure-state families only)</span>
            </label>
            <p id="qet-analysis-schmidt-hint" class="qet-generator__hint">
              This option is unavailable for mixed-state families.
            </p>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-schmidt" hidden></p>
          </div>
        </div>

        <div class="qet-generator__analysis-group" role="group" aria-labelledby="qet-analysis-certificates-title">
          <h3 id="qet-analysis-certificates-title">Classify separability</h3>
          <div class="qet-generator__checks">
            <label class="qet-generator__check">
              <input id="qet-analysis-pipeline" type="checkbox" checked>
              <span>Certificate-aware pipeline</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-pipeline" hidden></p>

            <label class="qet-generator__check">
              <input id="qet-analysis-separability" type="checkbox">
              <span>Structured separability report</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-separability" hidden></p>

            <div class="qet-generator__field qet-generator__nested-field">
              <label for="qet-separability-profile">Separability strategy profile</label>
              <select id="qet-separability-profile" disabled aria-describedby="qet-separability-profile-hint">
                <option value="ppt_only">PPT only — quickest necessary test</option>
                <option value="fast_detection">Fast entanglement detection</option>
                <option value="balanced_core" selected>Balanced deterministic core checks</option>
                <option value="separable_ball">Prioritize the sufficient separable-ball test</option>
              </select>
              <p id="qet-separability-profile-hint" class="qet-generator__hint">
                Profiles expand to a fixed, documented sequence of core
                strategies. They never load an optional solver.
              </p>
              <p class="qet-generator__field-error" data-qet-error-for="qet-separability-profile" hidden></p>
            </div>

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
          </div>
        </div>

        <div class="qet-generator__analysis-group" role="group" aria-labelledby="qet-analysis-reporting-title">
          <h3 id="qet-analysis-reporting-title">Measures and integration</h3>
          <div class="qet-generator__checks">
            <label class="qet-generator__check">
              <input id="qet-analysis-measures" type="checkbox" checked>
              <span>Purity, entropy, negativity, and logarithmic negativity</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-measures" hidden></p>

            <label class="qet-generator__check">
              <input id="qet-analysis-backend-status" type="checkbox">
              <span>Report optional-backend readiness without loading a backend</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-backend-status" hidden></p>
          </div>
        </div>

        <details class="qet-generator__advanced">
          <summary>Advanced tolerances</summary>
          <div class="qet-generator__field-grid qet-generator__tolerances">
            <div class="qet-generator__field">
              <label for="qet-analysis-atol">Absolute tolerance</label>
              <input id="qet-analysis-atol" type="text" inputmode="decimal" value="1e-12" aria-describedby="qet-tolerance-hint">
              <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-atol" hidden></p>
            </div>
            <div class="qet-generator__field">
              <label for="qet-analysis-rtol">Relative tolerance</label>
              <input id="qet-analysis-rtol" type="text" inputmode="decimal" value="0" aria-describedby="qet-tolerance-hint">
              <p class="qet-generator__field-error" data-qet-error-for="qet-analysis-rtol" hidden></p>
            </div>
          </div>
          <p id="qet-tolerance-hint" class="qet-generator__hint">
            Inputs are never normalized, clipped, symmetrized, or otherwise repaired.
          </p>
        </details>
      </fieldset>

      <fieldset class="qet-generator__section">
        <legend>4. Shape and save the example</legend>
        <div class="qet-generator__analysis-group" role="group" aria-labelledby="qet-output-options-title">
          <h3 id="qet-output-options-title">Generated output</h3>
          <div class="qet-generator__checks">
            <label class="qet-generator__check">
              <input id="qet-output-result-helpers" type="checkbox">
              <span>Use concise result helpers such as <code>conclusion</code> and <code>explain</code></span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-output-result-helpers" hidden></p>

            <label class="qet-generator__check">
              <input id="qet-output-include-assertions" type="checkbox">
              <span>Include safe structural assertions in the generated example</span>
            </label>
            <p class="qet-generator__field-error" data-qet-error-for="qet-output-include-assertions" hidden></p>
          </div>
        </div>

        <div class="qet-generator__analysis-group" role="group" aria-labelledby="qet-config-portability-title">
          <h3 id="qet-config-portability-title">Portable configuration</h3>
          <p id="qet-config-portability-hint" class="qet-generator__hint">
            Copy or download a canonical, versioned JSON configuration, or load
            one up to 16 KiB. Loading never evaluates code and does not generate
            Julia until you review the settings and press Generate.
          </p>
          <div class="qet-generator__config-actions">
            <button id="qet-generator-config-copy" class="qet-generator__button" type="button">
              Copy configuration
            </button>
            <button id="qet-generator-config-download" class="qet-generator__button" type="button">
              Download configuration
            </button>
          </div>
          <div class="qet-generator__field qet-generator__file-field">
            <label for="qet-generator-config-load">Load a JSON configuration</label>
            <input id="qet-generator-config-load" type="file" accept=".json,application/json" aria-describedby="qet-config-portability-hint">
            <p class="qet-generator__field-error" data-qet-error-for="qet-generator-config-load" hidden></p>
          </div>
          <p id="qet-generator-config-status" class="qet-generator__status qet-generator__config-status" aria-live="polite"></p>
        </div>
      </fieldset>

      <div id="qet-generator-errors" class="qet-generator__errors" role="alert" hidden></div>

      <div class="qet-generator__form-actions">
        <button class="qet-generator__button qet-generator__button--primary" type="submit">
          Generate Julia code
        </button>
        <button class="qet-generator__button" type="reset">Reset to default preset</button>
      </div>
    </form>

    <section class="qet-generator__output" aria-labelledby="qet-generator-output-title">
      <div class="qet-generator__output-heading">
        <div>
          <p class="qet-generator__eyebrow">Review and run</p>
          <h3 id="qet-generator-output-title">Generated Julia example</h3>
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
      <pre class="qet-generator__code" tabindex="0" aria-label="Generated Julia code" aria-describedby="qet-generator-status"><code id="qet-generator-code"># Loading generator assets…
</code></pre>
    </section>
  </div>
</div>
`;

    function byId(id) {
        return document.getElementById(id);
    }

    function setText(element, text) {
        if (element) {
            element.textContent = text;
        }
    }

    function appendOption(select, value, label) {
        var option = document.createElement("option");
        option.value = value;
        option.textContent = label;
        select.appendChild(option);
    }

    function populateSelectors() {
        var presetSelect = byId("qet-generator-preset");
        var familySelect = byId("qet-state-family");
        appendOption(presetSelect, "custom", "Custom selection");
        presetSelect.options[0].disabled = true;

        var presets = generator.presetList();
        for (var presetIndex = 0; presetIndex < presets.length; presetIndex += 1) {
            appendOption(
                presetSelect,
                presets[presetIndex].id,
                presets[presetIndex].label
            );
        }

        var families = generator.familyList();
        for (var familyIndex = 0; familyIndex < families.length; familyIndex += 1) {
            appendOption(
                familySelect,
                families[familyIndex].id,
                families[familyIndex].label
            );
        }
    }

    function fieldValue(id) {
        var element = byId(id);
        return element ? element.value : "";
    }

    function checked(id) {
        var element = byId(id);
        return Boolean(element && element.checked);
    }

    function readForm() {
        return {
            version: generator.schemaVersion,
            family: fieldValue("qet-state-family"),
            params: {
                dimA: fieldValue("qet-product-dim-a"),
                dimB: fieldValue("qet-product-dim-b"),
                indexA: fieldValue("qet-product-index-a"),
                indexB: fieldValue("qet-product-index-b"),
                bellIndex: fieldValue("qet-bell-index"),
                mixtureDimension: fieldValue("qet-mixture-dimension"),
                ghzDimension: fieldValue("qet-ghz-dimension"),
                ghzParties: fieldValue("qet-ghz-parties"),
                ghzCut: fieldValue("qet-ghz-cut"),
                dickeParties: fieldValue("qet-dicke-parties"),
                dickeExcitations: fieldValue("qet-dicke-excitations"),
                dickeCut: fieldValue("qet-dicke-cut"),
                isotropicDimension: fieldValue("qet-isotropic-dimension"),
                isotropicAlpha: fieldValue("qet-isotropic-alpha"),
                wernerDimension: fieldValue("qet-werner-dimension"),
                wernerAlpha: fieldValue("qet-werner-alpha"),
                horodeckiA: fieldValue("qet-horodecki-a"),
                horodeckiDims: fieldValue("qet-horodecki-dims"),
                symmetricP: fieldValue("qet-symmetric-p"),
                symmetricPhase: fieldValue("qet-symmetric-phase"),
                randomDimA: fieldValue("qet-random-dim-a"),
                randomDimB: fieldValue("qet-random-dim-b"),
                randomRank: fieldValue("qet-random-rank"),
                randomSeed: fieldValue("qet-random-seed"),
                randomDistribution: fieldValue("qet-random-distribution"),
                randomReal: checked("qet-random-real"),
            },
            analysis: {
                validation: checked("qet-analysis-validation"),
                marginals: checked("qet-analysis-marginals"),
                schmidt: checked("qet-analysis-schmidt"),
                pipeline: checked("qet-analysis-pipeline"),
                separability: checked("qet-analysis-separability"),
                ppt: checked("qet-analysis-ppt"),
                realignment: checked("qet-analysis-realignment"),
                reduction: checked("qet-analysis-reduction"),
                separable_ball: checked("qet-analysis-separable-ball"),
                measures: checked("qet-analysis-measures"),
                backend_status: checked("qet-analysis-backend-status"),
                ppt_witness: checked("qet-analysis-ppt-witness"),
            },
            options: {
                separabilityProfile: fieldValue("qet-separability-profile"),
            },
            output: {
                result_helpers: checked("qet-output-result-helpers"),
                include_assertions: checked("qet-output-include-assertions"),
            },
            atol: fieldValue("qet-analysis-atol"),
            rtol: fieldValue("qet-analysis-rtol"),
        };
    }

    function assignValue(id, value) {
        var element = byId(id);
        if (element && value !== undefined && value !== null) {
            element.value = String(value);
        }
    }

    function assignChecked(id, value) {
        var element = byId(id);
        if (element) {
            element.checked = Boolean(value);
        }
    }

    function applyConfig(config) {
        if (!config) {
            return;
        }
        applyingConfig = true;
        assignValue("qet-state-family", config.family);
        var params = config.params || {};
        assignValue("qet-product-dim-a", params.dimA);
        assignValue("qet-product-dim-b", params.dimB);
        assignValue("qet-product-index-a", params.indexA);
        assignValue("qet-product-index-b", params.indexB);
        assignValue("qet-bell-index", params.bellIndex);
        assignValue("qet-mixture-dimension", params.mixtureDimension);
        assignValue("qet-ghz-dimension", params.ghzDimension);
        assignValue("qet-ghz-parties", params.ghzParties);
        assignValue("qet-ghz-cut", params.ghzCut);
        assignValue("qet-dicke-parties", params.dickeParties);
        assignValue("qet-dicke-excitations", params.dickeExcitations);
        assignValue("qet-dicke-cut", params.dickeCut);
        assignValue("qet-isotropic-dimension", params.isotropicDimension);
        assignValue("qet-isotropic-alpha", params.isotropicAlpha);
        assignValue("qet-werner-dimension", params.wernerDimension);
        assignValue("qet-werner-alpha", params.wernerAlpha);
        assignValue("qet-horodecki-a", params.horodeckiA);
        assignValue("qet-horodecki-dims", params.horodeckiDims);
        assignValue("qet-symmetric-p", params.symmetricP);
        assignValue("qet-symmetric-phase", params.symmetricPhase);
        assignValue("qet-random-dim-a", params.randomDimA);
        assignValue("qet-random-dim-b", params.randomDimB);
        assignValue("qet-random-rank", params.randomRank);
        assignValue("qet-random-seed", params.randomSeed);
        assignValue("qet-random-distribution", params.randomDistribution);
        assignChecked("qet-random-real", params.randomReal);

        var analysis = config.analysis || {};
        assignChecked("qet-analysis-validation", analysis.validation);
        assignChecked("qet-analysis-marginals", analysis.marginals);
        assignChecked("qet-analysis-schmidt", analysis.schmidt);
        assignChecked("qet-analysis-pipeline", analysis.pipeline);
        assignChecked("qet-analysis-separability", analysis.separability);
        assignChecked("qet-analysis-ppt", analysis.ppt);
        assignChecked("qet-analysis-realignment", analysis.realignment);
        assignChecked("qet-analysis-reduction", analysis.reduction);
        assignChecked("qet-analysis-separable-ball", analysis.separable_ball);
        assignChecked("qet-analysis-measures", analysis.measures);
        assignChecked("qet-analysis-backend-status", analysis.backend_status);
        assignChecked("qet-analysis-ppt-witness", analysis.ppt_witness);
        var options = config.options || {};
        assignValue(
            "qet-separability-profile",
            options.separabilityProfile || "balanced_core"
        );
        var output = config.output || {};
        assignChecked("qet-output-result-helpers", output.result_helpers);
        assignChecked("qet-output-include-assertions", output.include_assertions);
        assignValue("qet-analysis-atol", config.atol);
        assignValue("qet-analysis-rtol", config.rtol);
        applyingConfig = false;
        updateConditionalFields();
    }

    function updateNumericBounds() {
        var dimA = Number(fieldValue("qet-product-dim-a"));
        var dimB = Number(fieldValue("qet-product-dim-b"));
        var ghzParties = Number(fieldValue("qet-ghz-parties"));
        var dickeParties = Number(fieldValue("qet-dicke-parties"));
        var isotropicDimension = Number(fieldValue("qet-isotropic-dimension"));
        var randomDimA = Number(fieldValue("qet-random-dim-a"));
        var randomDimB = Number(fieldValue("qet-random-dim-b"));

        if (isFinite(dimA)) {
            byId("qet-product-index-a").max = String(Math.max(0, dimA - 1));
        }
        if (isFinite(dimB)) {
            byId("qet-product-index-b").max = String(Math.max(0, dimB - 1));
        }
        if (isFinite(ghzParties)) {
            byId("qet-ghz-cut").max = String(Math.max(1, ghzParties - 1));
        }
        if (isFinite(dickeParties)) {
            byId("qet-dicke-cut").max = String(Math.max(1, dickeParties - 1));
            byId("qet-dicke-excitations").max = String(Math.max(0, dickeParties));
        }
        if (isFinite(isotropicDimension) && isotropicDimension >= 2) {
            byId("qet-isotropic-alpha").min = String(
                -1 / (isotropicDimension * isotropicDimension - 1)
            );
        }
        if (
            isFinite(randomDimA) &&
            isFinite(randomDimB) &&
            randomDimA >= 2 &&
            randomDimB >= 2
        ) {
            byId("qet-random-rank").max = String(randomDimA * randomDimB);
        }
    }

    function isPureFamily(family) {
        return (
            family === "product_basis" ||
            family === "bell" ||
            family === "ghz" ||
            family === "dicke"
        );
    }

    function positiveInteger(id) {
        var value = Number(fieldValue(id));
        if (!isFinite(value) || value <= 0 || Math.floor(value) !== value) {
            return null;
        }
        return value;
    }

    function formatBytes(bytes) {
        if (!isFinite(bytes) || bytes < 0) {
            return "unknown";
        }
        if (bytes < 1024) {
            return Math.ceil(bytes) + " B";
        }
        if (bytes < 1024 * 1024) {
            return (bytes / 1024).toFixed(bytes < 10 * 1024 ? 1 : 0) + " KiB";
        }
        return (bytes / (1024 * 1024)).toFixed(bytes < 10 * 1024 * 1024 ? 1 : 0) + " MiB";
    }

    function resourceDimensions(family) {
        var localDimension;
        var parties;
        var cut;
        var dimA;
        var dimB;

        if (family === "product_basis") {
            dimA = positiveInteger("qet-product-dim-a");
            dimB = positiveInteger("qet-product-dim-b");
        } else if (family === "bell") {
            dimA = 2;
            dimB = 2;
        } else if (family === "diagonal_mixture") {
            localDimension = positiveInteger("qet-mixture-dimension");
            dimA = localDimension;
            dimB = localDimension;
        } else if (family === "ghz") {
            localDimension = positiveInteger("qet-ghz-dimension");
            parties = positiveInteger("qet-ghz-parties");
            cut = positiveInteger("qet-ghz-cut");
            if (
                localDimension !== null &&
                parties !== null &&
                cut !== null &&
                cut < parties
            ) {
                dimA = Math.pow(localDimension, cut);
                dimB = Math.pow(localDimension, parties - cut);
            }
        } else if (family === "dicke") {
            parties = positiveInteger("qet-dicke-parties");
            cut = positiveInteger("qet-dicke-cut");
            if (parties !== null && cut !== null && cut < parties) {
                dimA = Math.pow(2, cut);
                dimB = Math.pow(2, parties - cut);
            }
        } else if (family === "isotropic") {
            localDimension = positiveInteger("qet-isotropic-dimension");
            dimA = localDimension;
            dimB = localDimension;
        } else if (family === "werner") {
            localDimension = positiveInteger("qet-werner-dimension");
            dimA = localDimension;
            dimB = localDimension;
        } else if (family === "horodecki") {
            if (fieldValue("qet-horodecki-dims") === "2x4") {
                dimA = 2;
                dimB = 4;
            } else {
                dimA = 3;
                dimB = 3;
            }
        } else if (family === "symmetric_sappt_ghz5") {
            dimA = 4;
            dimB = 8;
        } else if (family === "tiles_upb") {
            dimA = 3;
            dimB = 3;
        } else if (family === "random_density") {
            dimA = positiveInteger("qet-random-dim-a");
            dimB = positiveInteger("qet-random-dim-b");
        }

        if (dimA === null || dimB === null || dimA === undefined || dimB === undefined) {
            return null;
        }
        return [dimA, dimB];
    }

    function updateResourcePlan() {
        var plan = byId("qet-generator-resource-plan");
        if (!plan) {
            return;
        }
        var family = fieldValue("qet-state-family");
        var dimensions = resourceDimensions(family);
        plan.classList.remove("qet-generator__resource-plan--warning");
        if (!dimensions) {
            setText(
                plan,
                "Planning estimate unavailable until all dimensions and cut positions are positive integers."
            );
            return;
        }

        var totalDimension = dimensions[0] * dimensions[1];
        var denseEntries = totalDimension * totalDimension;
        var bytesPerEntry =
            family === "random_density" && checked("qet-random-real") ? 8 : 16;
        var description =
            (isPureFamily(family) ? "Pure-state family" : "Mixed-state family") +
            " · bipartition " +
            dimensions[0] +
            " × " +
            dimensions[1] +
            " · Hilbert dimension " +
            totalDimension +
            " · dense state " +
            denseEntries +
            " entries (about " +
            formatBytes(denseEntries * bytesPerEntry) +
            ").";

        var dimensionCap = family === "random_density" ? 64 : 256;
        if (totalDimension > dimensionCap) {
            description +=
                " This exceeds this family's documentation-generator dimension cap of " +
                dimensionCap +
                ".";
            plan.classList.add("qet-generator__resource-plan--warning");
        } else if (
            checked("qet-analysis-separability") &&
            totalDimension > 64
        ) {
            description +=
                " The structured separability report is capped at Hilbert dimension 64.";
            plan.classList.add("qet-generator__resource-plan--warning");
        } else {
            description +=
                " This is a planning estimate; generation performs the authoritative checks.";
        }
        setText(plan, description);
    }

    function updateConditionalFields() {
        var family = fieldValue("qet-state-family");
        var groups = root.querySelectorAll("[data-qet-family]");
        for (var index = 0; index < groups.length; index += 1) {
            var visible = groups[index].getAttribute("data-qet-family") === family;
            groups[index].hidden = !visible;
            var controls = groups[index].querySelectorAll("input, select");
            for (var controlIndex = 0; controlIndex < controls.length; controlIndex += 1) {
                controls[controlIndex].disabled = !visible;
            }
        }

        var pptSelected = checked("qet-analysis-ppt");
        var witness = byId("qet-analysis-ppt-witness");
        witness.disabled = !pptSelected;
        if (!pptSelected) {
            witness.checked = false;
        }

        var schmidt = byId("qet-analysis-schmidt");
        var mixedFamily = !isPureFamily(family);
        schmidt.disabled = mixedFamily;
        if (mixedFamily) {
            schmidt.checked = false;
        }

        var separabilitySelected = checked("qet-analysis-separability");
        byId("qet-separability-profile").disabled = !separabilitySelected;
        updateNumericBounds();
        updateResourcePlan();
    }

    function rememberBaseDescriptions() {
        var controls = root.querySelectorAll("input, select");
        for (var index = 0; index < controls.length; index += 1) {
            controls[index].setAttribute(
                "data-qet-base-describedby",
                controls[index].getAttribute("aria-describedby") || ""
            );
        }
    }

    function clearErrors() {
        var invalid = root.querySelectorAll("[aria-invalid='true']");
        for (var invalidIndex = 0; invalidIndex < invalid.length; invalidIndex += 1) {
            invalid[invalidIndex].removeAttribute("aria-invalid");
            var baseDescription = invalid[invalidIndex].getAttribute(
                "data-qet-base-describedby"
            );
            if (baseDescription) {
                invalid[invalidIndex].setAttribute(
                    "aria-describedby",
                    baseDescription
                );
            } else {
                invalid[invalidIndex].removeAttribute("aria-describedby");
            }
        }
        var localErrors = root.querySelectorAll("[data-qet-error-for]");
        for (var errorIndex = 0; errorIndex < localErrors.length; errorIndex += 1) {
            localErrors[errorIndex].hidden = true;
            localErrors[errorIndex].textContent = "";
        }
        var summary = byId("qet-generator-errors");
        while (summary.firstChild) {
            summary.removeChild(summary.firstChild);
        }
        summary.hidden = true;
    }

    function localFieldId(coreField) {
        var aliases = {
            "options-separability-profile": "qet-separability-profile",
            "analysis-separability-profile": "qet-separability-profile",
            "output-result-helpers": "qet-output-result-helpers",
            "output-include-assertions": "qet-output-include-assertions",
            "config": "qet-generator-config-load",
            "configuration": "qet-generator-config-load",
            "config-load": "qet-generator-config-load",
            "config-file": "qet-generator-config-load",
            "config-json": "qet-generator-config-load",
            "config-text": "qet-generator-config-load",
        };
        if (aliases[coreField]) {
            return aliases[coreField];
        }
        return "qet-" + coreField;
    }

    function focusControl(event) {
        var targetId = event.currentTarget.getAttribute("data-qet-focus");
        var control = byId(targetId);
        if (!control) {
            return;
        }
        event.preventDefault();
        var parent = control.parentNode;
        while (parent && parent !== root) {
            if (parent.tagName && parent.tagName.toLowerCase() === "details") {
                parent.open = true;
            }
            parent = parent.parentNode;
        }
        control.focus();
    }

    function describedByWithError(control, errorId) {
        var baseDescription =
            control.getAttribute("data-qet-base-describedby") || "";
        return baseDescription ? baseDescription + " " + errorId : errorId;
    }

    function showErrors(errors, focusFirst) {
        var summary = byId("qet-generator-errors");
        var heading = document.createElement("p");
        heading.textContent = "Please correct the following:";
        summary.appendChild(heading);
        var list = document.createElement("ul");
        var firstControl = null;

        for (var index = 0; index < errors.length; index += 1) {
            var item = document.createElement("li");
            var control = byId(localFieldId(errors[index].field));
            var localError = root.querySelector(
                "[data-qet-error-for='" + localFieldId(errors[index].field) + "']"
            );
            if (control) {
                var link = document.createElement("a");
                link.href = "#" + control.id;
                link.setAttribute("data-qet-focus", control.id);
                link.textContent = errors[index].message;
                link.addEventListener("click", focusControl);
                item.appendChild(link);
                control.setAttribute("aria-invalid", "true");
                if (localError) {
                    localError.id = localFieldId(errors[index].field) + "-error";
                    localError.textContent = errors[index].message;
                    localError.hidden = false;
                    control.setAttribute(
                        "aria-describedby",
                        describedByWithError(control, localError.id)
                    );
                }
                if (!firstControl) {
                    firstControl = control;
                }
            } else {
                item.textContent = errors[index].message;
            }
            list.appendChild(item);
        }
        summary.appendChild(list);
        summary.hidden = false;
        if (focusFirst && firstControl) {
            var parent = firstControl.parentNode;
            while (parent && parent !== root) {
                if (parent.tagName && parent.tagName.toLowerCase() === "details") {
                    parent.open = true;
                }
                parent = parent.parentNode;
            }
            firstControl.focus();
        }
    }

    function clearElement(element) {
        while (element.firstChild) {
            element.removeChild(element.firstChild);
        }
    }

    function summaryChip(text) {
        var chip = document.createElement("span");
        chip.className = "qet-generator__chip";
        chip.textContent = text;
        return chip;
    }

    function renderSummary(result) {
        var summary = byId("qet-generator-summary");
        clearElement(summary);
        summary.appendChild(summaryChip(result.summary.family));
        summary.appendChild(summaryChip(result.summary.stateKind + " state"));
        summary.appendChild(
            summaryChip(result.summary.dimensions.join(" × ") + " bipartition")
        );
        summary.appendChild(
            summaryChip("Hilbert dimension " + result.summary.totalDimension)
        );
        for (var index = 0; index < result.summary.analyses.length; index += 1) {
            summary.appendChild(summaryChip(result.summary.analyses[index]));
        }
    }

    function renderNotices(notices) {
        var container = byId("qet-generator-notices");
        clearElement(container);
        if (notices.length === 0) {
            container.hidden = true;
            return;
        }
        var list = document.createElement("ul");
        for (var index = 0; index < notices.length; index += 1) {
            var item = document.createElement("li");
            item.textContent = notices[index];
            list.appendChild(item);
        }
        container.appendChild(list);
        container.hidden = false;
    }

    function renderResult(result, focusErrors) {
        clearErrors();
        var copyButton = byId("qet-generator-copy");
        var downloadButton = byId("qet-generator-download");
        var output = root.querySelector(".qet-generator__output");
        var codeContainer = root.querySelector(".qet-generator__code");
        output.classList.remove("qet-generator__output--stale");
        codeContainer.setAttribute("aria-label", "Generated Julia code");
        if (!result.ok) {
            lastResult = null;
            showErrors(result.errors, focusErrors);
            setText(
                byId("qet-generator-code"),
                "# Fix the highlighted inputs, then generate again.\n"
            );
            clearElement(byId("qet-generator-summary"));
            renderNotices([]);
            copyButton.disabled = true;
            downloadButton.disabled = true;
            setText(
                byId("qet-generator-status"),
                "Code was not generated because some inputs are invalid."
            );
            return;
        }

        lastResult = result;
        setText(byId("qet-generator-code"), result.code);
        renderSummary(result);
        renderNotices(result.notices);
        copyButton.disabled = false;
        downloadButton.disabled = false;
        setText(
            byId("qet-generator-status"),
            "Generated " + result.filename + ". The code has not been executed in this page."
        );
    }

    function markOutputStale(message) {
        lastResult = null;
        byId("qet-generator-copy").disabled = true;
        byId("qet-generator-download").disabled = true;
        root
            .querySelector(".qet-generator__output")
            .classList.add("qet-generator__output--stale");
        root
            .querySelector(".qet-generator__code")
            .setAttribute(
                "aria-label",
                "Previously generated Julia code; stale settings require regeneration"
            );
        setText(
            byId("qet-generator-status"),
            message ||
                "Settings changed. Generate again before copying or downloading Julia code."
        );
    }

    function generateAndRender(focusErrors) {
        renderResult(generator.generate(readForm()), focusErrors);
    }

    function fallbackCopy(text) {
        var activeElement = document.activeElement;
        var selectionStart =
            activeElement && typeof activeElement.selectionStart === "number"
                ? activeElement.selectionStart
                : null;
        var selectionEnd =
            activeElement && typeof activeElement.selectionEnd === "number"
                ? activeElement.selectionEnd
                : null;
        var textarea = document.createElement("textarea");
        textarea.value = text;
        textarea.setAttribute("readonly", "");
        textarea.style.position = "fixed";
        textarea.style.opacity = "0";
        document.body.appendChild(textarea);
        textarea.select();
        var copied = false;
        try {
            copied = document.execCommand("copy");
        } finally {
            document.body.removeChild(textarea);
            if (activeElement && activeElement.focus) {
                activeElement.focus();
                if (
                    selectionStart !== null &&
                    selectionEnd !== null &&
                    activeElement.setSelectionRange
                ) {
                    activeElement.setSelectionRange(selectionStart, selectionEnd);
                }
            }
        }
        return copied;
    }

    function copyText(text, onSuccess, onFailure) {
        var fallback = function () {
            if (fallbackCopy(text)) {
                onSuccess();
            } else {
                onFailure();
            }
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(text).then(onSuccess, fallback);
        } else {
            fallback();
        }
    }

    function copyCode() {
        if (!lastResult) {
            return;
        }
        var onSuccess = function () {
            setText(byId("qet-generator-status"), "Generated Julia code copied.");
        };
        var onFailure = function () {
            setText(
                byId("qet-generator-status"),
                "Copy was blocked by the browser; select the code and copy it manually."
            );
        };
        copyText(lastResult.code, onSuccess, onFailure);
    }

    function downloadText(text, filename, mediaType) {
        var blob = new Blob([text], { type: mediaType });
        var url = URL.createObjectURL(blob);
        var link = document.createElement("a");
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        window.setTimeout(function () {
            URL.revokeObjectURL(url);
        }, 0);
    }

    function downloadCode() {
        if (!lastResult) {
            return;
        }
        downloadText(
            lastResult.code,
            lastResult.filename,
            "text/x-julia;charset=utf-8"
        );
        setText(
            byId("qet-generator-status"),
            "Downloaded " + lastResult.filename + "."
        );
    }

    function setConfigStatus(message) {
        setText(byId("qet-generator-config-status"), message);
    }

    function serializedFormConfig(focusErrors) {
        clearErrors();
        var serialized = generator.serializeConfig(readForm());
        if (!serialized.ok) {
            showErrors(serialized.errors, focusErrors);
            setConfigStatus(
                "Configuration was not exported because some settings are invalid."
            );
            return null;
        }
        return serialized;
    }

    function copyConfiguration() {
        var serialized = serializedFormConfig(true);
        if (!serialized) {
            return;
        }
        copyText(
            serialized.text,
            function () {
                setConfigStatus("Canonical JSON configuration copied.");
            },
            function () {
                setConfigStatus(
                    "Copy was blocked by the browser. Download the JSON configuration instead."
                );
            }
        );
    }

    function downloadConfiguration() {
        var serialized = serializedFormConfig(true);
        if (!serialized) {
            return;
        }
        downloadText(
            serialized.text,
            CONFIG_FILENAME,
            "application/json;charset=utf-8"
        );
        setConfigStatus("Downloaded " + CONFIG_FILENAME + ".");
    }

    function showConfigLoadError(message) {
        clearErrors();
        showErrors(
            [
                {
                    field: "generator-config-load",
                    message: message,
                },
            ],
            true
        );
        setConfigStatus(message);
    }

    function applyParsedConfiguration(parsed) {
        if (!parsed.ok) {
            clearErrors();
            var loadErrors = [];
            for (var index = 0; index < parsed.errors.length; index += 1) {
                loadErrors.push({
                    field: "generator-config-load",
                    message: parsed.errors[index].message,
                });
            }
            showErrors(loadErrors, true);
            setConfigStatus(
                "The JSON configuration was not loaded. Current settings were preserved."
            );
            return;
        }
        clearErrors();
        applyConfig(parsed.config);
        byId("qet-generator-preset").value = "custom";
        setText(
            byId("qet-preset-description"),
            "Loaded configuration. Review its bounded settings, then generate the Julia example."
        );
        updateConditionalFields();
        markOutputStale(
            "Configuration loaded. Review the settings and generate again before exporting Julia code."
        );
        setConfigStatus(
            "Configuration loaded locally. Nothing was executed, uploaded, or persisted."
        );
    }

    function loadConfiguration(event) {
        var input = event.target;
        var file = input.files && input.files.length > 0 ? input.files[0] : null;
        if (!file) {
            return;
        }
        if (file.size > CONFIG_MAX_BYTES) {
            input.value = "";
            showConfigLoadError("Configuration files must be at most 16 KiB.");
            return;
        }

        var reader = new FileReader();
        reader.onerror = function () {
            input.value = "";
            showConfigLoadError("The selected configuration file could not be read.");
        };
        reader.onload = function () {
            var text = typeof reader.result === "string" ? reader.result : "";
            input.value = "";
            if (text.length > CONFIG_MAX_BYTES) {
                showConfigLoadError("Configuration files must be at most 16 KiB.");
                return;
            }
            applyParsedConfiguration(generator.parseConfig(text));
        };
        reader.readAsText(file, "utf-8");
    }

    function choosePreset(id) {
        var config = generator.presetConfig(id);
        if (!config) {
            return;
        }
        applyConfig(config);
        byId("qet-generator-preset").value = id;
        var presets = generator.presetList();
        for (var index = 0; index < presets.length; index += 1) {
            if (presets[index].id === id) {
                setText(byId("qet-preset-description"), presets[index].description);
                break;
            }
        }
        setConfigStatus("");
        generateAndRender(false);
    }

    function markCustom() {
        if (applyingConfig) {
            return;
        }
        byId("qet-generator-preset").value = "custom";
        setText(
            byId("qet-preset-description"),
            "Custom settings stay in this page only unless you explicitly copy or download the configuration."
        );
    }

    function settingChanged(event) {
        if (
            event.target.id === "qet-generator-preset" ||
            event.target.id === "qet-generator-config-load" ||
            applyingConfig
        ) {
            return;
        }
        markCustom();
        clearErrors();
        setConfigStatus("");
        updateConditionalFields();
        markOutputStale();
    }

    function bindEvents() {
        byId("qet-generator-preset").addEventListener("change", function (event) {
            if (event.target.value !== "custom") {
                choosePreset(event.target.value);
            }
        });
        form.addEventListener("submit", function (event) {
            event.preventDefault();
            updateConditionalFields();
            generateAndRender(true);
        });
        form.addEventListener("reset", function (event) {
            event.preventDefault();
            choosePreset("separable_product");
            setText(
                byId("qet-generator-status"),
                "Reset to the certified pure-product preset."
            );
        });
        form.addEventListener("input", settingChanged);
        form.addEventListener("change", settingChanged);
        byId("qet-generator-copy").addEventListener("click", copyCode);
        byId("qet-generator-download").addEventListener("click", downloadCode);
        byId("qet-generator-config-copy").addEventListener(
            "click",
            copyConfiguration
        );
        byId("qet-generator-config-download").addEventListener(
            "click",
            downloadConfiguration
        );
        byId("qet-generator-config-load").addEventListener(
            "change",
            loadConfiguration
        );
    }

    function supportsGeneratorCore(candidate) {
        if (!candidate || candidate.schemaVersion !== 2) {
            return false;
        }
        var requiredMethods = [
            "normalizeConfig",
            "generate",
            "presetConfig",
            "presetList",
            "familyList",
            "serializeConfig",
            "parseConfig",
        ];
        for (var index = 0; index < requiredMethods.length; index += 1) {
            if (typeof candidate[requiredMethods[index]] !== "function") {
                return false;
            }
        }
        return true;
    }

    function showUnavailableFallback() {
        var controls = root.querySelectorAll("input, select, button");
        for (var index = 0; index < controls.length; index += 1) {
            controls[index].disabled = true;
        }
        root.classList.add("qet-generator--unavailable");
        setText(
            byId("qet-generator-status"),
            "The complete schema-2 generator assets did not load. Reload the page or use the static fallback below."
        );
        setText(
            byId("qet-generator-code"),
            "# Generator unavailable. Use the static fallback below.\n"
        );
        root
            .querySelector(".qet-generator__code")
            .setAttribute("aria-label", "Generator unavailable");
    }

    function initialize() {
        var mountHeading = byId("Interactive-generator");
        if (!mountHeading) {
            return;
        }
        mountHeading.insertAdjacentHTML("afterend", GENERATOR_MARKUP);
        root = byId("qet-code-generator");
        generator = window.QETCodeGenerator;
        form = byId("qet-generator-form");
        if (!form || !supportsGeneratorCore(generator)) {
            showUnavailableFallback();
            return;
        }

        root.setAttribute("data-generator-schema", String(generator.schemaVersion));
        populateSelectors();
        rememberBaseDescriptions();
        bindEvents();
        choosePreset("separable_product");
        root.classList.add("qet-generator--ready");
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initialize);
    } else {
        initialize();
    }
})();

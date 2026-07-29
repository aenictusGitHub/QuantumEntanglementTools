(function () {
    "use strict";

    var generator = null;
    var root = null;
    var form = null;
    var lastResult = null;
    var applyingConfig = false;

    // Keep the full widget markup out of the Markdown source: GitHub does not
    // interpret Documenter's `@raw html` blocks, while built docs load this asset.
    var GENERATOR_MARKUP = `
<div id="qet-code-generator" class="qet-generator" data-generator-schema="1">
  <div class="qet-generator__local-note" role="note">
    <strong>Local and deterministic.</strong>
    This page generates text in your browser. It does not run Julia, send the
    selected parameters to a server, or save them after the page is reloaded.
  </div>


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
            },
            analysis: {
                pipeline: checked("qet-analysis-pipeline"),
                ppt: checked("qet-analysis-ppt"),
                realignment: checked("qet-analysis-realignment"),
                reduction: checked("qet-analysis-reduction"),
                separable_ball: checked("qet-analysis-separable-ball"),
                measures: checked("qet-analysis-measures"),
                ppt_witness: checked("qet-analysis-ppt-witness"),
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

        var analysis = config.analysis || {};
        assignChecked("qet-analysis-pipeline", analysis.pipeline);
        assignChecked("qet-analysis-ppt", analysis.ppt);
        assignChecked("qet-analysis-realignment", analysis.realignment);
        assignChecked("qet-analysis-reduction", analysis.reduction);
        assignChecked("qet-analysis-separable-ball", analysis.separable_ball);
        assignChecked("qet-analysis-measures", analysis.measures);
        assignChecked("qet-analysis-ppt-witness", analysis.ppt_witness);
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
        updateNumericBounds();
    }

    function clearErrors() {
        var invalid = root.querySelectorAll("[aria-invalid='true']");
        for (var invalidIndex = 0; invalidIndex < invalid.length; invalidIndex += 1) {
            invalid[invalidIndex].removeAttribute("aria-invalid");
            invalid[invalidIndex].removeAttribute("aria-describedby");
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
        return "qet-" + coreField;
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
            item.textContent = errors[index].message;
            list.appendChild(item);

            var control = byId(localFieldId(errors[index].field));
            var localError = root.querySelector(
                "[data-qet-error-for='" + localFieldId(errors[index].field) + "']"
            );
            if (control) {
                control.setAttribute("aria-invalid", "true");
                if (localError) {
                    localError.id = localFieldId(errors[index].field) + "-error";
                    localError.textContent = errors[index].message;
                    localError.hidden = false;
                    control.setAttribute("aria-describedby", localError.id);
                }
                if (!firstControl) {
                    firstControl = control;
                }
            }
        }
        summary.appendChild(list);
        summary.hidden = false;
        if (focusFirst && firstControl) {
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

    function generateAndRender(focusErrors) {
        renderResult(generator.generate(readForm()), focusErrors);
    }

    function fallbackCopy(text) {
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
        }
        return copied;
    }

    function copyCode() {
        if (!lastResult) {
            return;
        }
        var onSuccess = function () {
            setText(byId("qet-generator-status"), "Generated Julia code copied.");
        };
        var onFailure = function () {
            if (fallbackCopy(lastResult.code)) {
                onSuccess();
            } else {
                setText(
                    byId("qet-generator-status"),
                    "Copy was blocked by the browser; select the code and copy it manually."
                );
            }
        };

        if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(lastResult.code).then(onSuccess, onFailure);
        } else {
            onFailure();
        }
    }

    function downloadCode() {
        if (!lastResult) {
            return;
        }
        var blob = new Blob([lastResult.code], { type: "text/x-julia;charset=utf-8" });
        var url = URL.createObjectURL(blob);
        var link = document.createElement("a");
        link.href = url;
        link.download = lastResult.filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        window.setTimeout(function () {
            URL.revokeObjectURL(url);
        }, 0);
        setText(
            byId("qet-generator-status"),
            "Downloaded " + lastResult.filename + "."
        );
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
        generateAndRender(false);
    }

    function markCustom() {
        if (applyingConfig) {
            return;
        }
        byId("qet-generator-preset").value = "custom";
        setText(
            byId("qet-preset-description"),
            "Custom settings are kept only in this page until it is reloaded."
        );
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
        form.addEventListener("input", function (event) {
            if (event.target.id === "qet-generator-preset") {
                return;
            }
            markCustom();
            clearErrors();
            updateConditionalFields();
        });
        form.addEventListener("change", function (event) {
            if (event.target.id === "qet-generator-preset") {
                return;
            }
            markCustom();
            clearErrors();
            updateConditionalFields();
        });
        byId("qet-generator-copy").addEventListener("click", copyCode);
        byId("qet-generator-download").addEventListener("click", downloadCode);
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
        if (!generator || !form) {
            setText(
                byId("qet-generator-status"),
                "The generator assets did not load. Reload the page or use the static examples below."
            );
            return;
        }

        populateSelectors();
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

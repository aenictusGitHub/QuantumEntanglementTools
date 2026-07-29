(function () {
    "use strict";

    var generator = null;
    var root = null;
    var form = null;
    var lastResult = null;
    var applyingConfig = false;

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
        root = byId("qet-code-generator");
        if (!root) {
            return;
        }
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

(function (root, factory) {
    "use strict";

    var api = factory();
    if (typeof module === "object" && module.exports) {
        module.exports = api;
    } else {
        root.QETCodeGeneratorCases = api;
    }
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
    "use strict";

    var REQUIRED_DOM_IDS = [
        "qet-code-generator",
        "qet-generator-form",
        "qet-generator-preset",
        "qet-state-family",
        "qet-analysis-pipeline",
        "qet-analysis-ppt",
        "qet-analysis-ppt-witness",
        "qet-analysis-realignment",
        "qet-analysis-reduction",
        "qet-analysis-separable-ball",
        "qet-analysis-measures",
        "qet-analysis-atol",
        "qet-analysis-rtol",
        "qet-generator-errors",
        "qet-generator-summary",
        "qet-generator-notices",
        "qet-generator-status",
        "qet-generator-code",
        "qet-generator-copy",
        "qet-generator-download",
    ];

    function clone(value) {
        return JSON.parse(JSON.stringify(value));
    }

    function baseConfig(family, params) {
        return {
            version: 1,
            family: family,
            params: params,
            analysis: {
                pipeline: true,
                ppt: false,
                realignment: false,
                reduction: false,
                separable_ball: false,
                measures: false,
                ppt_witness: false,
            },
            atol: 1e-12,
            rtol: 0,
        };
    }

    function hasErrorCode(result, code) {
        for (var index = 0; index < result.errors.length; index += 1) {
            if (result.errors[index].code === code) {
                return true;
            }
        }
        return false;
    }

    function moduleText(name, generatedCode, assertions) {
        return (
            "module " +
            name +
            "\n\n" +
            generatedCode +
            "\n" +
            assertions.join("\n") +
            "\n\nend # module " +
            name +
            "\n"
        );
    }

    function run(core) {
        var failures = [];
        var checks = 0;

        function check(condition, label) {
            checks += 1;
            if (!condition) {
                failures.push(label);
            }
        }

        var families = [
            baseConfig("product_basis", { dimA: 2, dimB: 3, indexA: 1, indexB: 2 }),
            baseConfig("bell", { bellIndex: 3 }),
            baseConfig("diagonal_mixture", { mixtureDimension: 3 }),
            baseConfig("ghz", { ghzDimension: 2, ghzParties: 3, ghzCut: 1 }),
            baseConfig("dicke", { dickeParties: 4, dickeExcitations: 2, dickeCut: 2 }),
            baseConfig("isotropic", { isotropicDimension: 3, isotropicAlpha: 0.5 }),
            baseConfig("werner", { wernerDimension: 3, wernerAlpha: 0.8 }),
            baseConfig("horodecki", { horodeckiA: 0.3, horodeckiDims: "3x3" }),
            baseConfig("symmetric_sappt_ghz5", {
                symmetricP: 121 / 125,
                symmetricPhase: 0,
            }),
        ];

        for (var familyIndex = 0; familyIndex < families.length; familyIndex += 1) {
            var input = families[familyIndex];
            var before = JSON.stringify(input);
            var first = core.generate(input);
            var second = core.generate(input);
            check(first.ok, "family " + input.family + " must generate");
            check(first.code === second.code, "family " + input.family + " must be deterministic");
            check(
                JSON.stringify(input) === before,
                "family " + input.family + " must not mutate its input"
            );
            check(
                first.code.charAt(first.code.length - 1) === "\n",
                "family " + input.family + " must end with a newline"
            );
            check(
                !/(undefined|NaN|Infinity)/.test(first.code),
                "family " + input.family + " must not emit non-finite placeholders"
            );
        }

        var presetIds = core.presetList().map(function (preset) {
            return preset.id;
        });
        check(presetIds.length === 5, "five curated presets must be available");
        for (var presetIndex = 0; presetIndex < presetIds.length; presetIndex += 1) {
            var presetConfig = core.presetConfig(presetIds[presetIndex]);
            var presetResult = core.generate(presetConfig);
            check(presetResult.ok, "preset " + presetIds[presetIndex] + " must generate");
            check(
                core.presetConfig(presetIds[presetIndex]) !== presetConfig,
                "preset " + presetIds[presetIndex] + " must return a defensive copy"
            );
        }
        check(core.presetConfig("missing") === null, "unknown presets must be rejected");

        var malicious = baseConfig("bell", { bellIndex: "0; error(\"injected\")" });
        var maliciousResult = core.generate(malicious);
        check(!maliciousResult.ok, "Julia injection text must be rejected");
        check(
            hasErrorCode(maliciousResult, "invalid_number"),
            "Julia injection must produce invalid_number"
        );

        var whitespace = baseConfig("bell", { bellIndex: "   " });
        check(!core.generate(whitespace).ok, "whitespace-only numbers must be rejected");

        var empty = baseConfig("bell", { bellIndex: "" });
        check(!core.generate(empty).ok, "empty numeric fields must be rejected");

        var fractional = baseConfig("bell", { bellIndex: 1.5 });
        check(
            hasErrorCode(core.generate(fractional), "integer_required"),
            "fractional integer fields must be rejected"
        );

        var unsupported = baseConfig("bell", { bellIndex: 0 });
        unsupported.version = 2;
        check(
            hasErrorCode(core.generate(unsupported), "unsupported_schema"),
            "unsupported schemas must be rejected"
        );

        var outOfRangeIsotropic = baseConfig("isotropic", {
            isotropicDimension: 3,
            isotropicAlpha: -0.2,
        });
        check(
            hasErrorCode(core.generate(outOfRangeIsotropic), "number_below_minimum"),
            "the dimension-dependent isotropic positivity range must be enforced"
        );

        var oversized = baseConfig("ghz", {
            ghzDimension: 4,
            ghzParties: 8,
            ghzCut: 4,
        });
        check(
            hasErrorCode(core.generate(oversized), "resource_limit"),
            "dense GHZ output must obey the resource cap"
        );

        var noAnalysis = baseConfig("bell", { bellIndex: 0 });
        noAnalysis.analysis.pipeline = false;
        check(
            hasErrorCode(core.generate(noAnalysis), "analysis_required"),
            "at least one analysis must be selected"
        );

        var orphanWitness = baseConfig("bell", { bellIndex: 0 });
        orphanWitness.analysis.ppt_witness = true;
        check(
            hasErrorCode(core.generate(orphanWitness), "witness_requires_ppt"),
            "a PPT witness must require the PPT route"
        );

        var sparseFamilies = [
            ["ghz", "sparse_output=false"],
            ["dicke", "sparse_output=false"],
            ["isotropic", "sparse_output=false"],
            ["werner", "sparse_output=false"],
        ];
        for (var sparseIndex = 0; sparseIndex < sparseFamilies.length; sparseIndex += 1) {
            var matching = null;
            for (var searchIndex = 0; searchIndex < families.length; searchIndex += 1) {
                if (families[searchIndex].family === sparseFamilies[sparseIndex][0]) {
                    matching = families[searchIndex];
                    break;
                }
            }
            var denseCode = core.generate(matching).code;
            check(
                denseCode.indexOf(sparseFamilies[sparseIndex][1]) !== -1,
                sparseFamilies[sparseIndex][0] + " must explicitly request dense output"
            );
        }

        var modules = [];

        var product = core.generate(core.presetConfig("separable_product"));
        modules.push(
            moduleText("GeneratedProduct", product.code, [
                "@assert report.status === :separable",
                "@assert report.certified",
                "@assert report.certificate_kind === :pure_product_decomposition",
            ])
        );

        var bell = core.generate(core.presetConfig("bell_witness"));
        modules.push(
            moduleText("GeneratedBell", bell.code, [
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :pure_state_schmidt_rank",
                "@assert ppt_result.status === CriterionEntanglementDetected",
                "@assert ppt_witness !== nothing",
                "@assert isapprox(state_negativity, 0.5; atol=1e-12, rtol=0)",
                "@assert isapprox(log_negativity, 1; atol=1e-12, rtol=0)",
            ])
        );

        var mixture = core.generate(core.presetConfig("separable_mixture"));
        modules.push(
            moduleText("GeneratedMixture", mixture.code, [
                "@assert report.status === :unknown",
                "@assert !report.certified",
                "@assert ball_result.status === :separable_certified",
            ])
        );

        var horodecki = core.generate(core.presetConfig("horodecki_ppt_entangled"));
        modules.push(
            moduleText("GeneratedHorodecki", horodecki.code, [
                "@assert ppt_result.status !== CriterionEntanglementDetected",
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :realignment_cross_norm_violation",
            ])
        );

        var ghzConfig = baseConfig("ghz", {
            ghzDimension: 2,
            ghzParties: 3,
            ghzCut: 1,
        });
        var ghz = core.generate(ghzConfig);
        modules.push(
            moduleText("GeneratedGHZ", ghz.code, [
                "@assert dims == (2, 4)",
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :pure_state_schmidt_rank",
            ])
        );

        var dickeConfig = baseConfig("dicke", {
            dickeParties: 4,
            dickeExcitations: 2,
            dickeCut: 2,
        });
        var dicke = core.generate(dickeConfig);
        modules.push(
            moduleText("GeneratedDicke", dicke.code, [
                "@assert dims == (4, 4)",
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :pure_state_schmidt_rank",
            ])
        );

        var isotropicConfig = baseConfig("isotropic", {
            isotropicDimension: 3,
            isotropicAlpha: 0.5,
        });
        var isotropic = core.generate(isotropicConfig);
        modules.push(
            moduleText("GeneratedIsotropic", isotropic.code, [
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :negative_partial_transpose_witness",
            ])
        );

        var wernerConfig = baseConfig("werner", {
            wernerDimension: 3,
            wernerAlpha: 0.8,
        });
        var werner = core.generate(wernerConfig);
        modules.push(
            moduleText("GeneratedWerner", werner.code, [
                "@assert report.status === :entangled",
                "@assert report.certified",
                "@assert report.certificate_kind === :negative_partial_transpose_witness",
            ])
        );

        var symmetric = core.generate(core.presetConfig("symmetric_sappt_witness"));
        modules.push(
            moduleText("GeneratedSymmetricSAPPT", symmetric.code, [
                "@assert isapprox(real(tr(rho)), 1; atol=1e-12, rtol=0)",
                "@assert w5_expectation < 0",
            ])
        );

        var bundle =
            "# Generated code-generator smoke bundle. Do not edit.\n\n" +
            modules.join("\n") +
            '\nprintln("Generated code smoke examples passed.")\n';

        return {
            checks: checks,
            failures: failures,
            bundle: bundle,
            requiredDomIds: REQUIRED_DOM_IDS.slice(),
        };
    }

    return { run: run, requiredDomIds: REQUIRED_DOM_IDS.slice() };
});

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
        "qet-analysis-validation",
        "qet-analysis-marginals",
        "qet-analysis-schmidt",
        "qet-analysis-pipeline",
        "qet-analysis-separability",
        "qet-analysis-ppt",
        "qet-analysis-ppt-witness",
        "qet-analysis-realignment",
        "qet-analysis-reduction",
        "qet-analysis-separable-ball",
        "qet-analysis-measures",
        "qet-analysis-backend-status",
        "qet-analysis-atol",
        "qet-analysis-rtol",
        "qet-separability-profile",
        "qet-output-result-helpers",
        "qet-output-include-assertions",
        "qet-random-dim-a",
        "qet-random-dim-b",
        "qet-random-rank",
        "qet-random-seed",
        "qet-random-distribution",
        "qet-random-real",
        "qet-generator-errors",
        "qet-generator-summary",
        "qet-generator-notices",
        "qet-generator-resource-plan",
        "qet-generator-status",
        "qet-generator-code",
        "qet-generator-copy",
        "qet-generator-download",
        "qet-generator-config-copy",
        "qet-generator-config-download",
        "qet-generator-config-load",
        "qet-generator-config-status",
    ];

    function clone(value) {
        return JSON.parse(JSON.stringify(value));
    }

    function baseConfig(family, params) {
        return {
            version: 2,
            family: family,
            params: params || {},
            analysis: {
                pipeline: true,
                ppt: false,
                realignment: false,
                reduction: false,
                separable_ball: false,
                measures: false,
                validation: false,
                marginals: false,
                schmidt: false,
                separability: false,
                backend_status: false,
                ppt_witness: false,
            },
            atol: 1e-12,
            rtol: 0,
            options: {
                separabilityProfile: "balanced_core",
            },
            output: {
                result_helpers: false,
                include_assertions: false,
            },
        };
    }

    function legacyConfig() {
        return {
            version: 1,
            family: "bell",
            params: { bellIndex: 0 },
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
        if (!result || !result.errors) {
            return false;
        }
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

    function validateDocumentation(pageSource, uiSource) {
        var failures = [];
        var templateMatch = /var GENERATOR_MARKUP = `([\s\S]*?)`;/m.exec(uiSource);
        if (!templateMatch) {
            failures.push("UI adapter must define the generator markup template");
        }

        var renderedSource = pageSource + "\n" + (templateMatch ? templateMatch[1] : "");
        if (pageSource.indexOf("## Interactive generator") !== -1) {
            renderedSource += '\n<h2 id="Interactive-generator"></h2>';
        }
        for (var idIndex = 0; idIndex < REQUIRED_DOM_IDS.length; idIndex += 1) {
            var id = REQUIRED_DOM_IDS[idIndex];
            if (renderedSource.indexOf('id="' + id + '"') === -1) {
                failures.push("rendered generator is missing required DOM id " + id);
            }
        }

        var referencedIdPattern = /byId\("([^"]+)"\)/g;
        var referencedIdMatch;
        while ((referencedIdMatch = referencedIdPattern.exec(uiSource)) !== null) {
            if (renderedSource.indexOf('id="' + referencedIdMatch[1] + '"') === -1) {
                failures.push(
                    "UI adapter references missing DOM id " + referencedIdMatch[1]
                );
            }
        }

        if (pageSource.indexOf("```@raw html") !== -1) {
            failures.push("documentation page must not expose raw HTML in GitHub previews");
        }
        if (pageSource.indexOf("## Interactive generator") === -1) {
            failures.push("documentation page must provide the generator mount heading");
        }
        if (pageSource.indexOf('id="qet-generator-form"') !== -1) {
            failures.push("interactive form markup must live in the UI asset");
        }
        if (uiSource.indexOf('byId("Interactive-generator")') === -1) {
            failures.push("UI adapter must target the generated mount-heading id");
        }
        if (pageSource.indexOf("Viewing this source on GitHub?") === -1) {
            failures.push("documentation page must explain the GitHub source preview");
        }
        if (pageSource.indexOf("/actions/workflows/docs.yml") === -1) {
            failures.push("documentation page must link to the Documentation workflow");
        }
        if (
            pageSource.indexOf(
                "`documentation-<run-id>-<run-attempt>` artifact"
            ) === -1
        ) {
            failures.push(
                "documentation page must identify the run-specific docs artifact"
            );
        }
        if (pageSource.indexOf("`code_generator/index.html`") === -1) {
            failures.push("documentation page must identify the artifact entry point");
        }
        if (
            pageSource.indexOf(
                "julia --startup-file=no scripts/build_docs.jl"
            ) === -1
        ) {
            failures.push("documentation page must provide the local docs build command");
        }
        if (pageSource.indexOf("[executable tutorials](tutorials.md)") === -1) {
            failures.push("documentation page must link to executable tutorials");
        }
        if (
            pageSource.indexOf("## License of generated artifacts") === -1 ||
            pageSource.indexOf("`SPDX-License-Identifier`") === -1 ||
            pageSource.indexOf("descriptive metadata") === -1
        ) {
            failures.push(
                "documentation page must explain generated-code and configuration licensing"
            );
        }

        return failures;
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

        check(core.schemaVersion === 2, "the public generator schema must be version 2");
        check(
            typeof core.serializeConfig === "function" &&
                typeof core.parseConfig === "function",
            "configuration portability helpers must be public"
        );
        check(core.familyList().length === 11, "all eleven state families must be listed");

        var families = [
            baseConfig("product_basis", {
                dimA: 2,
                dimB: 3,
                indexA: 1,
                indexB: 2,
            }),
            baseConfig("bell", { bellIndex: 3 }),
            baseConfig("diagonal_mixture", { mixtureDimension: 3 }),
            baseConfig("ghz", { ghzDimension: 2, ghzParties: 3, ghzCut: 1 }),
            baseConfig("dicke", {
                dickeParties: 4,
                dickeExcitations: 2,
                dickeCut: 2,
            }),
            baseConfig("isotropic", {
                isotropicDimension: 3,
                isotropicAlpha: 0.5,
            }),
            baseConfig("werner", { wernerDimension: 3, wernerAlpha: 0.8 }),
            baseConfig("horodecki", {
                horodeckiA: 0.3,
                horodeckiDims: "3x3",
            }),
            baseConfig("symmetric_sappt_ghz5", {
                symmetricP: 121 / 125,
                symmetricPhase: 0,
            }),
            baseConfig("tiles_upb", {}),
            baseConfig("random_density", {
                randomDimA: 2,
                randomDimB: 3,
                randomRank: 3,
                randomSeed: 20260731,
                randomDistribution: "hilbert_schmidt",
                randomReal: false,
            }),
        ];

        for (var familyIndex = 0; familyIndex < families.length; familyIndex += 1) {
            var input = families[familyIndex];
            var before = JSON.stringify(input);
            var first = core.generate(input);
            var second = core.generate(input);
            var hasTemplateAttribution =
                first.code.indexOf("# SPDX-FileCopyrightText: 2026 John Martin") !==
                    -1 &&
                first.code.indexOf("# SPDX-License-Identifier: BSD-3-Clause") !==
                    -1 &&
                first.code.indexOf(
                    "QuantumEntanglementTools/blob/main/LICENSE"
                ) !== -1;
            var hasPaperAttribution =
                input.family !== "symmetric_sappt_ghz5" ||
                (first.code.indexOf("10.1103/PhysRevA.111.042418") !== -1 &&
                    first.code.indexOf("does not copy or rerun the source SDP") !==
                        -1);
            check(
                first.ok && hasTemplateAttribution && hasPaperAttribution,
                "family " +
                    input.family +
                    " must generate with required template and source attribution"
            );
            check(
                first.code === second.code,
                "family " + input.family + " must be deterministic"
            );
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

        var requiredPresetIds = [
            "separable_product",
            "bell_witness",
            "separable_mixture",
            "horodecki_ppt_entangled",
            "symmetric_sappt_witness",
            "tiles_bound_entanglement",
            "seeded_random_diagnostics",
            "ghz_schmidt",
            "separability_strategy_comparison",
        ];
        var presetIds = core.presetList().map(function (preset) {
            return preset.id;
        });
        check(presetIds.length >= 9, "at least nine curated presets must be available");
        for (
            var requiredPresetIndex = 0;
            requiredPresetIndex < requiredPresetIds.length;
            requiredPresetIndex += 1
        ) {
            check(
                presetIds.indexOf(requiredPresetIds[requiredPresetIndex]) !== -1,
                "required preset " +
                    requiredPresetIds[requiredPresetIndex] +
                    " must remain available"
            );
        }
        for (var presetIndex = 0; presetIndex < presetIds.length; presetIndex += 1) {
            var preset = core.presetConfig(presetIds[presetIndex]);
            var presetResult = core.generate(preset);
            check(presetResult.ok, "preset " + presetIds[presetIndex] + " must generate");
            check(
                core.presetConfig(presetIds[presetIndex]) !== preset,
                "preset " + presetIds[presetIndex] + " must return a defensive copy"
            );
        }
        check(core.presetConfig("missing") === null, "unknown presets must be rejected");

        var legacy = legacyConfig();
        var legacyBefore = JSON.stringify(legacy);
        var migrated = core.normalizeConfig(legacy);
        check(migrated.ok, "schema-v1 configurations must migrate");
        check(migrated.config.version === 2, "migration must produce schema v2");
        check(
            !migrated.config.analysis.validation &&
                !migrated.config.analysis.marginals &&
                !migrated.config.analysis.schmidt &&
                !migrated.config.analysis.separability &&
                !migrated.config.analysis.backend_status,
            "migration must leave all new analysis switches disabled"
        );
        check(
            migrated.config.output.result_helpers === false &&
                migrated.config.output.include_assertions === false,
            "migration must leave both new output switches disabled"
        );
        check(
            migrated.config.options.separabilityProfile === "balanced_core",
            "migration must select the dependency-free balanced profile"
        );
        check(
            JSON.stringify(legacy) === legacyBefore,
            "migration must not mutate schema-v1 input"
        );
        var missingVersion = legacyConfig();
        delete missingVersion.version;
        check(
            core.parseConfig(JSON.stringify(missingVersion)).ok,
            "a missing schema version must follow the v1 migration"
        );
        var future = legacyConfig();
        future.version = 3;
        check(
            hasErrorCode(core.parseConfig(JSON.stringify(future)), "unsupported_schema"),
            "future schemas must be rejected"
        );

        var portableSource = core.presetConfig("ghz_schmidt");
        var portableBefore = JSON.stringify(portableSource);
        var serialized = core.serializeConfig(portableSource);
        check(serialized.ok, "valid configurations must serialize");
        check(
            JSON.stringify(portableSource) === portableBefore,
            "serialization must not mutate input"
        );
        check(
            serialized.text.indexOf('"derived"') === -1,
            "portable JSON must omit derived fields"
        );
        check(
            serialized.text.charAt(serialized.text.length - 1) === "\n",
            "canonical JSON must end with one newline"
        );
        var parsed = core.parseConfig(serialized.text);
        check(parsed.ok, "canonical JSON must parse");
        check(
            core.generate(parsed.config).code === core.generate(portableSource).code,
            "serialize/parse round trips must preserve generated code"
        );
        check(
            core.serializeConfig(parsed.config).text === serialized.text,
            "serialized JSON must be canonical"
        );
        check(
            hasErrorCode(core.parseConfig("{"), "invalid_json"),
            "malformed JSON must be rejected"
        );
        check(
            hasErrorCode(core.parseConfig("[]"), "invalid_json"),
            "non-object JSON must be rejected"
        );
        check(
            hasErrorCode(
                core.parseConfig('{"version":2,"family":"bell","intruder":true}'),
                "unknown_config_key"
            ),
            "unknown top-level keys must be rejected"
        );
        var nestedIntruder = core.serializeConfig(baseConfig("bell", { bellIndex: 0 }));
        var nestedObject = JSON.parse(nestedIntruder.text);
        nestedObject.params.constructorName = "error";
        check(
            hasErrorCode(
                core.parseConfig(JSON.stringify(nestedObject)),
                "unknown_config_key"
            ),
            "unknown nested keys must be rejected"
        );
        var crossFamilyParameter = JSON.parse(nestedIntruder.text);
        crossFamilyParameter.params.randomSeed = 7;
        check(
            hasErrorCode(
                core.parseConfig(JSON.stringify(crossFamilyParameter)),
                "unknown_config_key"
            ),
            "parameters from another state family must not be silently discarded"
        );
        check(
            hasErrorCode(core.parseConfig(new Array(16386).join(" ")), "config_too_large"),
            "configuration text larger than 16 KiB must be rejected"
        );
        var multibyteOversized =
            '{"padding":"' + new Array(9001).join("\u00e9") + '"}';
        check(
            multibyteOversized.length < 16384 &&
                hasErrorCode(
                    core.parseConfig(multibyteOversized),
                    "config_too_large"
                ),
            "the 16 KiB limit must count UTF-8 bytes rather than UTF-16 code units"
        );

        var malicious = baseConfig("bell", {
            bellIndex: '0; error("injected")',
        });
        check(
            hasErrorCode(core.generate(malicious), "invalid_number"),
            "Julia injection text must be rejected as an invalid number"
        );
        check(
            !core.generate(baseConfig("bell", { bellIndex: "   " })).ok,
            "whitespace-only numbers must be rejected"
        );
        check(
            !core.generate(baseConfig("bell", { bellIndex: "" })).ok,
            "empty numeric fields must be rejected"
        );
        check(
            hasErrorCode(
                core.generate(baseConfig("bell", { bellIndex: 1.5 })),
                "integer_required"
            ),
            "fractional integer fields must be rejected"
        );
        check(
            hasErrorCode(
                core.generate(baseConfig("bell", { bellIndex: NaN })),
                "invalid_number"
            ) &&
                hasErrorCode(
                    core.generate(baseConfig("bell", { bellIndex: Infinity })),
                    "invalid_number"
                ),
            "non-finite JavaScript numbers must be rejected before emission"
        );
        check(
            hasErrorCode(
                core.generate(
                    baseConfig("isotropic", {
                        isotropicDimension: 3,
                        isotropicAlpha: -0.2,
                    })
                ),
                "number_below_minimum"
            ),
            "the dimension-dependent isotropic positivity range must be enforced"
        );
        check(
            hasErrorCode(
                core.generate(
                    baseConfig("ghz", {
                        ghzDimension: 4,
                        ghzParties: 8,
                        ghzCut: 4,
                    })
                ),
                "resource_limit"
            ),
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
        var mixedSchmidt = baseConfig("tiles_upb", {});
        mixedSchmidt.analysis.pipeline = false;
        mixedSchmidt.analysis.schmidt = true;
        check(
            hasErrorCode(core.generate(mixedSchmidt), "schmidt_requires_pure"),
            "Schmidt analysis must reject mixed-state families"
        );
        var largeSeparability = baseConfig("ghz", {
            ghzDimension: 2,
            ghzParties: 7,
            ghzCut: 3,
        });
        largeSeparability.analysis.pipeline = false;
        largeSeparability.analysis.separability = true;
        check(
            hasErrorCode(
                core.generate(largeSeparability),
                "separability_resource_limit"
            ),
            "composite separability analysis must obey its dimension cap"
        );

        var invalidRandomCases = [
            [
                {
                    randomDimA: 8,
                    randomDimB: 8,
                    randomRank: 65,
                    randomSeed: 0,
                    randomDistribution: "hilbert_schmidt",
                    randomReal: false,
                },
                "number_above_maximum",
            ],
            [
                {
                    randomDimA: 2,
                    randomDimB: 3,
                    randomRank: 3,
                    randomSeed: 4294967296,
                    randomDistribution: "hilbert_schmidt",
                    randomReal: false,
                },
                "number_above_maximum",
            ],
            [
                {
                    randomDimA: 2,
                    randomDimB: 3,
                    randomRank: 3,
                    randomSeed: 0,
                    randomDistribution: "wishart",
                    randomReal: false,
                },
                "invalid_choice",
            ],
            [
                {
                    randomDimA: 2,
                    randomDimB: 3,
                    randomRank: 3,
                    randomSeed: 0,
                    randomDistribution: "bures",
                    randomReal: "yes",
                },
                "boolean_required",
            ],
        ];
        for (
            var invalidRandomIndex = 0;
            invalidRandomIndex < invalidRandomCases.length;
            invalidRandomIndex += 1
        ) {
            check(
                hasErrorCode(
                    core.generate(
                        baseConfig(
                            "random_density",
                            invalidRandomCases[invalidRandomIndex][0]
                        )
                    ),
                    invalidRandomCases[invalidRandomIndex][1]
                ),
                "invalid random-state case " + invalidRandomIndex + " must fail"
            );
        }
        var randomDimensionLimit = baseConfig("random_density", {
            randomDimA: 8,
            randomDimB: 8,
            randomRank: 8,
            randomSeed: 0,
            randomDistribution: "bures",
            randomReal: true,
        });
        check(
            core.generate(randomDimensionLimit).ok,
            "the maximum permitted random-state dimension must generate"
        );

        var invalidProfile = baseConfig("bell", { bellIndex: 0 });
        invalidProfile.options.separabilityProfile = "full";
        check(
            hasErrorCode(core.generate(invalidProfile), "invalid_choice"),
            "unbounded or optional separability profiles must be rejected"
        );

        var profiles = {
            ppt_only: "(:ppt,)",
            fast_detection: "(:ppt, :realignment, :reduction)",
            balanced_core:
                "(:ppt, :realignment, :reduction, :separable_ball)",
            separable_ball: "(:separable_ball,)",
        };
        for (var profile in profiles) {
            if (Object.prototype.hasOwnProperty.call(profiles, profile)) {
                var profileConfig = baseConfig("diagonal_mixture", {
                    mixtureDimension: 2,
                });
                profileConfig.analysis.pipeline = false;
                profileConfig.analysis.separability = true;
                profileConfig.options.separabilityProfile = profile;
                var profileResult = core.generate(profileConfig);
                check(profileResult.ok, "profile " + profile + " must generate");
                check(
                    profileResult.code.indexOf(
                        "separability_strategies = " + profiles[profile]
                    ) !== -1,
                    "profile " + profile + " must emit its explicit tuple"
                );
                check(
                    profileResult.code.indexOf(
                        "available_separability_strategies()"
                    ) !== -1 &&
                        profileResult.code.indexOf("describe_strategy(strategy)") !==
                            -1,
                    "profile " + profile + " must emit API-drift checks"
                );
                check(
                    profileResult.code.indexOf(
                        "max_dense_entries=65_536, max_work=10_000_000"
                    ) !== -1,
                    "profile " + profile + " must emit fixed resource caps"
                );
                check(
                    profileResult.code.indexOf(
                        "if iszero(imag(tr(rho))) && real(tr(rho)) == 1"
                    ) !== -1,
                    "profile " + profile + " must guard exact real unit trace"
                );
                check(
                    profileResult.code.indexOf("not normalized or repaired") !== -1,
                    "profile " + profile + " must explain conservative skipping"
                );
            }
        }

        var helperConfig = baseConfig("bell", { bellIndex: 0 });
        helperConfig.analysis.ppt = true;
        helperConfig.output.result_helpers = true;
        helperConfig.output.include_assertions = true;
        var helperResult = core.generate(helperConfig);
        check(helperResult.ok, "result-helper mode must generate");
        check(
            helperResult.code.indexOf('": raw status="') !== -1 &&
                helperResult.code.indexOf("conclusion(result)") !== -1 &&
                helperResult.code.indexOf("is_conclusive(result)") !== -1 &&
                helperResult.code.indexOf("is_certified(result)") !== -1 &&
                helperResult.code.indexOf("explain(result)") !== -1,
            "helper mode must show native status and all conservative helpers"
        );
        check(
            helperResult.code.indexOf("@assert report.status") === -1 &&
                helperResult.code.indexOf("@assert conclusion(") === -1 &&
                helperResult.code.indexOf("@assert witness_expectation") === -1,
            "generated assertions must not claim parameter-dependent conclusions"
        );
        check(
            helperResult.code.indexOf("@assert size(rho)") !== -1 &&
                helperResult.code.indexOf("@assert report isa EntanglementReport") !==
                    -1,
            "assertion mode must emit structural checks"
        );

        var tilesCode = core.generate(
            core.presetConfig("tiles_bound_entanglement")
        ).code;
        check(
            tilesCode.indexOf("BigInt[") !== -1 &&
                tilesCode.indexOf("Rational{BigInt}") !== -1 &&
                tilesCode.indexOf("rho = Float64.(rho_exact)") !== -1,
            "Tiles must be constructed exactly before the Float64 copy"
        );
        check(
            tilesCode.indexOf(
                "tiles_complement = Matrix{Rational{BigInt}}(I, 9, 9) - tiles_projector"
            ) !== -1 &&
                tilesCode.indexOf(
                    "tiles_partial_transpose_exact = partial_transpose"
                ) !== -1,
            "Tiles must emit the exact complement and exact PPT calculation"
        );
        check(
            tilesCode.split("error isa InterruptException && rethrow()").length - 1 ===
                2,
            "generated separability and measure wrappers must preserve interrupts"
        );

        var randomCode = core.generate(
            core.presetConfig("seeded_random_diagnostics")
        ).code;
        check(
            randomCode.indexOf("using Random") !== -1 &&
                randomCode.indexOf("Xoshiro(random_seed)") !== -1,
            "random states must use an explicit local Xoshiro"
        );
        check(
            randomCode.indexOf("backend_readiness.entanglement.native.loaded") !==
                -1 &&
                randomCode.indexOf(
                    "@assert !backend_readiness.optimization.configured"
                ) !== -1,
            "backend diagnostics must expose only stable structural readiness checks"
        );

        var validationCode = core.generate(
            core.presetConfig("seeded_random_diagnostics")
        ).code;
        check(
            validationCode.indexOf("validation_report.valid") !== -1 &&
                validationCode.indexOf("@assert validation_report.valid") !== -1,
            "validation must report and structurally assert the valid field"
        );
        var ghzSchmidtCode = core.generate(
            core.presetConfig("ghz_schmidt")
        ).code;
        check(
            ghzSchmidtCode.indexOf("schmidt_decomposition(psi, dims") !== -1 &&
                ghzSchmidtCode.indexOf("schmidt_rank(psi, dims") !== -1,
            "the pure-state Schmidt route must emit both decomposition and rank"
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
            check(
                core.generate(matching).code.indexOf(sparseFamilies[sparseIndex][1]) !==
                    -1,
                sparseFamilies[sparseIndex][0] +
                    " must explicitly request dense output"
            );
        }

        var modules = [];
        var smokePresets = [
            ["GeneratedProduct", "separable_product", ["@assert size(rho) == (4, 4)"]],
            ["GeneratedBell", "bell_witness", ["@assert size(rho) == (4, 4)"]],
            [
                "GeneratedMixture",
                "separable_mixture",
                ["@assert size(rho) == (9, 9)"],
            ],
            [
                "GeneratedHorodecki",
                "horodecki_ppt_entangled",
                ["@assert size(rho) == (9, 9)"],
            ],
            [
                "GeneratedSymmetricSAPPT",
                "symmetric_sappt_witness",
                ["@assert size(rho) == (32, 32)"],
            ],
            [
                "GeneratedTiles",
                "tiles_bound_entanglement",
                [
                    "@assert tr(rho_exact) == 1",
                    "@assert tiles_partial_transpose_exact == rho_exact",
                ],
            ],
            [
                "GeneratedRandom",
                "seeded_random_diagnostics",
                [
                    "@assert size(rho) == (6, 6)",
                    "@assert state_purity !== nothing",
                    "@assert entropy_bits !== nothing",
                    "function generated_interrupt_probe()",
                    "    try",
                    '        evaluate_without_repair("Interrupt probe", () -> throw(InterruptException()))',
                    "    catch error",
                    "        return error isa InterruptException",
                    "    end",
                    "    return false",
                    "end",
                    "@assert generated_interrupt_probe()",
                ],
            ],
            [
                "GeneratedGHZSchmidt",
                "ghz_schmidt",
                ["@assert length(schmidt_result.coefficients) == 2"],
            ],
            [
                "GeneratedSeparabilityComparison",
                "separability_strategy_comparison",
                [
                    "@assert separability_strategies == (:ppt, :realignment, :reduction, :separable_ball)",
                ],
            ],
        ];
        for (
            var smokePresetIndex = 0;
            smokePresetIndex < smokePresets.length;
            smokePresetIndex += 1
        ) {
            var smokePreset = smokePresets[smokePresetIndex];
            var smokeResult = core.generate(core.presetConfig(smokePreset[1]));
            check(smokeResult.ok, "smoke preset " + smokePreset[1] + " must generate");
            modules.push(
                moduleText(smokePreset[0], smokeResult.code, smokePreset[2])
            );
        }

        var randomBures = baseConfig("random_density", {
            randomDimA: 2,
            randomDimB: 2,
            randomRank: 2,
            randomSeed: 0,
            randomDistribution: "bures",
            randomReal: true,
        });
        randomBures.analysis.pipeline = false;
        randomBures.analysis.validation = true;
        randomBures.output.include_assertions = true;

        var extraFamilies = [
            [
                "GeneratedDicke",
                baseConfig("dicke", {
                    dickeParties: 4,
                    dickeExcitations: 2,
                    dickeCut: 2,
                }),
                "@assert dims == (4, 4)",
            ],
            [
                "GeneratedIsotropic",
                baseConfig("isotropic", {
                    isotropicDimension: 3,
                    isotropicAlpha: 0.5,
                }),
                "@assert dims == (3, 3)",
            ],
            [
                "GeneratedWerner",
                baseConfig("werner", {
                    wernerDimension: 3,
                    wernerAlpha: 0.8,
                }),
                "@assert dims == (3, 3)",
            ],
            [
                "GeneratedRandomBuresReal",
                randomBures,
                "@assert size(rho) == (4, 4)",
            ],
        ];
        for (
            var extraFamilyIndex = 0;
            extraFamilyIndex < extraFamilies.length;
            extraFamilyIndex += 1
        ) {
            var extra = extraFamilies[extraFamilyIndex];
            var extraResult = core.generate(extra[1]);
            check(extraResult.ok, "extra smoke family " + extra[0] + " must generate");
            modules.push(moduleText(extra[0], extraResult.code, [extra[2]]));
        }

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

    return {
        run: run,
        requiredDomIds: REQUIRED_DOM_IDS.slice(),
        validateDocumentation: validateDocumentation,
    };
});

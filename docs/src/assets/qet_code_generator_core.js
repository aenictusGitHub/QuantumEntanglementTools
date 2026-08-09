(function (root, factory) {
    "use strict";

    var api = factory();
    if (typeof module === "object" && module.exports) {
        module.exports = api;
    } else {
        root.QETCodeGenerator = api;
    }
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
    "use strict";

    var SCHEMA_VERSION = 2;
    var MAX_CONFIG_BYTES = 16 * 1024;
    var FAMILY_ORDER = [
        "product_basis",
        "bell",
        "diagonal_mixture",
        "ghz",
        "dicke",
        "isotropic",
        "werner",
        "horodecki",
        "symmetric_sappt_ghz5",
        "tiles_upb",
        "random_density",
    ];
    var FAMILY_LABELS = {
        product_basis: "Computational-basis product state",
        bell: "Bell state",
        diagonal_mixture: "Explicit separable diagonal mixture",
        ghz: "GHZ state across a contiguous cut",
        dicke: "Dicke state across a contiguous cut",
        isotropic: "Isotropic state",
        werner: "Werner state",
        horodecki: "Horodecki state",
        symmetric_sappt_ghz5: "Five-qubit symmetric SAPPT family",
        tiles_upb: "Tiles-UPB complementary bound-entangled state",
        random_density: "Seeded random density matrix",
    };
    var ANALYSIS_ORDER = [
        "pipeline",
        "ppt",
        "realignment",
        "reduction",
        "separable_ball",
        "measures",
        "validation",
        "marginals",
        "schmidt",
        "separability",
        "backend_status",
    ];
    var ANALYSIS_LABELS = {
        pipeline: "certificate-aware pipeline",
        ppt: "PPT criterion",
        realignment: "realignment / CCNR",
        reduction: "reduction criterion",
        separable_ball: "separable-ball certificate",
        measures: "entropy and negativities",
        validation: "density-matrix validation",
        marginals: "reduced-state diagnostics",
        schmidt: "Schmidt decomposition",
        separability: "bounded separability strategy search",
        backend_status: "backend readiness",
    };
    var SEPARABILITY_PROFILES = {
        ppt_only: ["ppt"],
        fast_detection: ["ppt", "realignment", "reduction"],
        balanced_core: ["ppt", "realignment", "reduction", "separable_ball"],
        separable_ball: ["separable_ball"],
    };
    var FAMILY_PARAM_KEYS = {
        product_basis: ["dimA", "dimB", "indexA", "indexB"],
        bell: ["bellIndex"],
        diagonal_mixture: ["mixtureDimension"],
        ghz: ["ghzDimension", "ghzParties", "ghzCut"],
        dicke: ["dickeParties", "dickeExcitations", "dickeCut"],
        isotropic: ["isotropicDimension", "isotropicAlpha"],
        werner: ["wernerDimension", "wernerAlpha"],
        horodecki: ["horodeckiA", "horodeckiDims"],
        symmetric_sappt_ghz5: ["symmetricP", "symmetricPhase"],
        tiles_upb: [],
        random_density: [
            "randomDimA",
            "randomDimB",
            "randomRank",
            "randomSeed",
            "randomDistribution",
            "randomReal",
        ],
    };

    function clone(value) {
        if (Array.isArray(value)) {
            return value.map(function (item) {
                return clone(item);
            });
        }
        if (value && typeof value === "object") {
            var copied = {};
            for (var key in value) {
                if (hasOwn(value, key)) {
                    copied[key] = clone(value[key]);
                }
            }
            return copied;
        }
        return value;
    }

    function utf8ByteLength(text) {
        var bytes = 0;
        for (var index = 0; index < text.length; index += 1) {
            var code = text.charCodeAt(index);
            if (code <= 0x7f) {
                bytes += 1;
            } else if (code <= 0x7ff) {
                bytes += 2;
            } else if (
                code >= 0xd800 &&
                code <= 0xdbff &&
                index + 1 < text.length &&
                text.charCodeAt(index + 1) >= 0xdc00 &&
                text.charCodeAt(index + 1) <= 0xdfff
            ) {
                bytes += 4;
                index += 1;
            } else {
                bytes += 3;
            }
        }
        return bytes;
    }

    function hasOwn(object, key) {
        return Object.prototype.hasOwnProperty.call(object, key);
    }

    function objectOrEmpty(value) {
        return value && typeof value === "object" && !Array.isArray(value) ? value : {};
    }

    function rawConfig(family, params, analysis, atol, rtol, options, output) {
        return {
            version: SCHEMA_VERSION,
            family: family,
            params: params || {},
            analysis: analysis || {},
            atol: atol === undefined ? 1e-12 : atol,
            rtol: rtol === undefined ? 0 : rtol,
            options: options || { separabilityProfile: "balanced_core" },
            output: output || { result_helpers: false, include_assertions: false },
        };
    }

    function presetOutput() {
        return { result_helpers: true, include_assertions: true };
    }

    var PRESETS = {
        separable_product: {
            label: "Certified pure product",
            description: "Build |0> tensor |1> and obtain an exact Schmidt product certificate.",
            config: rawConfig(
                "product_basis",
                { dimA: 2, dimB: 2, indexA: 0, indexB: 1 },
                {
                    pipeline: true,
                    ppt: false,
                    realignment: false,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    ppt_witness: false,
                },
                undefined,
                undefined,
                undefined,
                presetOutput()
            ),
        },
        bell_witness: {
            label: "Bell state and PPT witness",
            description: "Detect Bell-state entanglement and construct a decomposable witness.",
            config: rawConfig(
                "bell",
                { bellIndex: 0 },
                {
                    pipeline: true,
                    ppt: true,
                    realignment: true,
                    reduction: true,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    ppt_witness: true,
                },
                undefined,
                undefined,
                undefined,
                presetOutput()
            ),
        },
        separable_mixture: {
            label: "Higher-dimensional separable mixture",
            description:
                "Compare inconclusive necessary tests with a sufficient separable-ball certificate.",
            config: rawConfig(
                "diagonal_mixture",
                { mixtureDimension: 3 },
                {
                    pipeline: true,
                    ppt: true,
                    realignment: true,
                    reduction: true,
                    separable_ball: true,
                    measures: true,
                    validation: true,
                    ppt_witness: false,
                },
                undefined,
                undefined,
                undefined,
                presetOutput()
            ),
        },
        horodecki_ppt_entangled: {
            label: "Horodecki PPT-entangled state",
            description:
                "Show that PPT can be inconclusive while realignment certifies entanglement.",
            config: rawConfig(
                "horodecki",
                { horodeckiA: 0.3, horodeckiDims: "3x3" },
                {
                    pipeline: true,
                    ppt: true,
                    realignment: true,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    ppt_witness: false,
                },
                1e-12,
                1e-10,
                undefined,
                presetOutput()
            ),
        },
        symmetric_sappt_witness: {
            label: "Five-qubit SAPPT witness",
            description:
                "Construct the paper's symmetric family and its rounded phase-matched witness.",
            config: rawConfig(
                "symmetric_sappt_ghz5",
                { symmetricP: 121 / 125, symmetricPhase: 0 },
                {
                    pipeline: true,
                    ppt: true,
                    realignment: true,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    ppt_witness: false,
                },
                1e-12,
                0,
                undefined,
                presetOutput()
            ),
        },
        tiles_bound_entanglement: {
            label: "Tiles bound entanglement",
            description:
                "Build the Tiles complement exactly, verify its density-matrix structure, and run bounded native detection.",
            config: rawConfig(
                "tiles_upb",
                {},
                {
                    pipeline: true,
                    ppt: true,
                    realignment: true,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    marginals: true,
                    schmidt: false,
                    separability: true,
                    backend_status: false,
                    ppt_witness: false,
                },
                1e-12,
                0,
                { separabilityProfile: "fast_detection" },
                presetOutput()
            ),
        },
        seeded_random_diagnostics: {
            label: "Seeded random-state diagnostics",
            description:
                "Generate a reproducible full-rank 2 x 3 density matrix and inspect validation, marginals, measures, and backend readiness.",
            config: rawConfig(
                "random_density",
                {
                    randomDimA: 2,
                    randomDimB: 3,
                    randomRank: 6,
                    randomSeed: 20260731,
                    randomDistribution: "hilbert_schmidt",
                    randomReal: false,
                },
                {
                    pipeline: false,
                    ppt: false,
                    realignment: false,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    marginals: true,
                    schmidt: false,
                    separability: false,
                    backend_status: true,
                    ppt_witness: false,
                },
                1e-12,
                1e-10,
                undefined,
                presetOutput()
            ),
        },
        ghz_schmidt: {
            label: "GHZ Schmidt decomposition",
            description:
                "Construct a three-qubit GHZ state and inspect its bipartite Schmidt data and marginals.",
            config: rawConfig(
                "ghz",
                { ghzDimension: 2, ghzParties: 3, ghzCut: 1 },
                {
                    pipeline: true,
                    ppt: false,
                    realignment: false,
                    reduction: false,
                    separable_ball: false,
                    measures: true,
                    validation: true,
                    marginals: true,
                    schmidt: true,
                    separability: false,
                    backend_status: false,
                    ppt_witness: false,
                },
                1e-12,
                0,
                undefined,
                presetOutput()
            ),
        },
        separability_strategy_comparison: {
            label: "Separability strategy comparison",
            description:
                "Compare one-sided criteria with a bounded, core-only separability certificate search.",
            config: rawConfig(
                "diagonal_mixture",
                { mixtureDimension: 3 },
                {
                    pipeline: false,
                    ppt: true,
                    realignment: true,
                    reduction: true,
                    separable_ball: true,
                    measures: false,
                    validation: true,
                    marginals: false,
                    schmidt: false,
                    separability: true,
                    backend_status: false,
                    ppt_witness: false,
                },
                1e-12,
                0,
                { separabilityProfile: "balanced_core" },
                presetOutput()
            ),
        },
    };

    function addError(errors, code, field, message) {
        errors.push({ code: code, field: field, message: message });
    }

    function finiteNumber(value, fallback, field, label, minimum, maximum, errors) {
        var candidate = value === undefined || value === null ? fallback : value;
        var number;
        if (typeof candidate === "number") {
            number = candidate;
        } else if (
            typeof candidate === "string" &&
            /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/.test(candidate.trim())
        ) {
            number = Number(candidate);
        } else {
            number = NaN;
        }
        if (typeof number !== "number" || !isFinite(number)) {
            addError(errors, "invalid_number", field, label + " must be a finite number.");
            return fallback;
        }
        if (minimum !== null && number < minimum) {
            addError(
                errors,
                "number_below_minimum",
                field,
                label + " must be at least " + minimum + "."
            );
        }
        if (maximum !== null && number > maximum) {
            addError(
                errors,
                "number_above_maximum",
                field,
                label + " must be at most " + maximum + "."
            );
        }
        return number;
    }

    function finiteInteger(value, fallback, field, label, minimum, maximum, errors) {
        var number = finiteNumber(value, fallback, field, label, minimum, maximum, errors);
        if (Math.floor(number) !== number) {
            addError(errors, "integer_required", field, label + " must be an integer.");
        }
        return number;
    }

    function enumValue(value, fallback, allowed, field, label, errors) {
        var candidate = value === undefined || value === null || value === "" ? fallback : value;
        if (typeof candidate !== "string" || allowed.indexOf(candidate) === -1) {
            addError(
                errors,
                "invalid_choice",
                field,
                label + " must be one of: " + allowed.join(", ") + "."
            );
            return fallback;
        }
        return candidate;
    }

    function booleanSetting(source, key, fallback, field, label, errors) {
        if (!hasOwn(source, key)) {
            return fallback;
        }
        if (typeof source[key] !== "boolean") {
            addError(
                errors,
                "boolean_required",
                field,
                label + " must be selected or cleared explicitly."
            );
            return fallback;
        }
        return source[key];
    }

    function booleanValue(source, key, fallback, errors) {
        return booleanSetting(
            source,
            key,
            fallback,
            "analysis-" + key.replace(/_/g, "-"),
            ANALYSIS_LABELS[key] || key,
            errors
        );
    }

    function migrateConfig(raw) {
        var source = clone(objectOrEmpty(raw));
        var version =
            !hasOwn(source, "version") || source.version === undefined
                ? 1
                : source.version;
        if (version !== 1 && version !== SCHEMA_VERSION) {
            return {
                ok: false,
                errors: [
                    {
                        code: "unsupported_schema",
                        field: "state-family",
                        message:
                            "This configuration uses an unsupported generator schema version.",
                    },
                ],
                config: source,
            };
        }
        if (version === SCHEMA_VERSION) {
            return { ok: true, errors: [], config: source };
        }

        var migratedAnalysis = clone(objectOrEmpty(source.analysis));
        migratedAnalysis.validation = false;
        migratedAnalysis.marginals = false;
        migratedAnalysis.schmidt = false;
        migratedAnalysis.separability = false;
        migratedAnalysis.backend_status = false;
        source.version = SCHEMA_VERSION;
        source.analysis = migratedAnalysis;
        source.options = { separabilityProfile: "balanced_core" };
        source.output = {
            result_helpers: false,
            include_assertions: false,
        };
        return { ok: true, errors: [], config: source };
    }

    function normalizeConfig(raw) {
        var migration = migrateConfig(raw);
        var errors = migration.errors.slice();
        var source = objectOrEmpty(migration.config);
        var family = source.family;
        if (typeof family !== "string" || FAMILY_ORDER.indexOf(family) === -1) {
            addError(
                errors,
                "unknown_family",
                "state-family",
                "Choose one of the supported state families."
            );
            family = "bell";
        }

        var rawParams = objectOrEmpty(source.params);
        var params = {};
        var totalDimension = 4;
        var dimensions = [2, 2];
        var stateKind = "pure";

        if (family === "product_basis") {
            params.dimA = finiteInteger(
                rawParams.dimA,
                2,
                "product-dim-a",
                "Subsystem A dimension",
                2,
                8,
                errors
            );
            params.dimB = finiteInteger(
                rawParams.dimB,
                2,
                "product-dim-b",
                "Subsystem B dimension",
                2,
                8,
                errors
            );
            params.indexA = finiteInteger(
                rawParams.indexA,
                0,
                "product-index-a",
                "Subsystem A basis label",
                0,
                Math.max(0, params.dimA - 1),
                errors
            );
            params.indexB = finiteInteger(
                rawParams.indexB,
                1,
                "product-index-b",
                "Subsystem B basis label",
                0,
                Math.max(0, params.dimB - 1),
                errors
            );
            totalDimension = params.dimA * params.dimB;
            dimensions = [params.dimA, params.dimB];
        } else if (family === "bell") {
            params.bellIndex = finiteInteger(
                rawParams.bellIndex,
                0,
                "bell-index",
                "Bell-state index",
                0,
                3,
                errors
            );
        } else if (family === "diagonal_mixture") {
            params.mixtureDimension = finiteInteger(
                rawParams.mixtureDimension,
                3,
                "mixture-dimension",
                "Local dimension",
                2,
                3,
                errors
            );
            if (params.mixtureDimension !== 2 && params.mixtureDimension !== 3) {
                addError(
                    errors,
                    "invalid_mixture_dimension",
                    "mixture-dimension",
                    "The explicit mixture is available in local dimension 2 or 3."
                );
            }
            totalDimension = params.mixtureDimension * params.mixtureDimension;
            dimensions = [params.mixtureDimension, params.mixtureDimension];
            stateKind = "mixed";
        } else if (family === "ghz") {
            params.ghzDimension = finiteInteger(
                rawParams.ghzDimension,
                2,
                "ghz-dimension",
                "GHZ local dimension",
                2,
                4,
                errors
            );
            params.ghzParties = finiteInteger(
                rawParams.ghzParties,
                3,
                "ghz-parties",
                "Number of parties",
                2,
                8,
                errors
            );
            params.ghzCut = finiteInteger(
                rawParams.ghzCut,
                1,
                "ghz-cut",
                "Parties on the left of the cut",
                1,
                Math.max(1, params.ghzParties - 1),
                errors
            );
            totalDimension = Math.pow(params.ghzDimension, params.ghzParties);
            dimensions = [
                Math.pow(params.ghzDimension, params.ghzCut),
                Math.pow(params.ghzDimension, params.ghzParties - params.ghzCut),
            ];
            if (totalDimension > 256) {
                addError(
                    errors,
                    "resource_limit",
                    "ghz-parties",
                    "The generated dense GHZ analysis is capped at total dimension 256."
                );
            }
        } else if (family === "dicke") {
            params.dickeParties = finiteInteger(
                rawParams.dickeParties,
                4,
                "dicke-parties",
                "Number of qubits",
                2,
                8,
                errors
            );
            params.dickeExcitations = finiteInteger(
                rawParams.dickeExcitations,
                2,
                "dicke-excitations",
                "Number of excitations",
                0,
                Math.max(0, params.dickeParties),
                errors
            );
            params.dickeCut = finiteInteger(
                rawParams.dickeCut,
                2,
                "dicke-cut",
                "Qubits on the left of the cut",
                1,
                Math.max(1, params.dickeParties - 1),
                errors
            );
            totalDimension = Math.pow(2, params.dickeParties);
            dimensions = [
                Math.pow(2, params.dickeCut),
                Math.pow(2, params.dickeParties - params.dickeCut),
            ];
        } else if (family === "isotropic") {
            params.isotropicDimension = finiteInteger(
                rawParams.isotropicDimension,
                3,
                "isotropic-dimension",
                "Isotropic local dimension",
                2,
                8,
                errors
            );
            var isotropicLower =
                -1 / (params.isotropicDimension * params.isotropicDimension - 1);
            params.isotropicAlpha = finiteNumber(
                rawParams.isotropicAlpha,
                0.5,
                "isotropic-alpha",
                "Isotropic alpha",
                isotropicLower,
                1,
                errors
            );
            totalDimension = params.isotropicDimension * params.isotropicDimension;
            dimensions = [params.isotropicDimension, params.isotropicDimension];
            stateKind = "mixed";
        } else if (family === "werner") {
            params.wernerDimension = finiteInteger(
                rawParams.wernerDimension,
                3,
                "werner-dimension",
                "Werner local dimension",
                2,
                8,
                errors
            );
            params.wernerAlpha = finiteNumber(
                rawParams.wernerAlpha,
                0.8,
                "werner-alpha",
                "Werner alpha",
                -1,
                1,
                errors
            );
            totalDimension = params.wernerDimension * params.wernerDimension;
            dimensions = [params.wernerDimension, params.wernerDimension];
            stateKind = "mixed";
        } else if (family === "horodecki") {
            params.horodeckiA = finiteNumber(
                rawParams.horodeckiA,
                0.3,
                "horodecki-a",
                "Horodecki a",
                0,
                1,
                errors
            );
            params.horodeckiDims = enumValue(
                rawParams.horodeckiDims,
                "3x3",
                ["3x3", "2x4"],
                "horodecki-dims",
                "Horodecki dimensions",
                errors
            );
            dimensions = params.horodeckiDims === "3x3" ? [3, 3] : [2, 4];
            totalDimension = dimensions[0] * dimensions[1];
            stateKind = "mixed";
        } else if (family === "symmetric_sappt_ghz5") {
            params.symmetricP = finiteNumber(
                rawParams.symmetricP,
                121 / 125,
                "symmetric-p",
                "Symmetric-family mixing weight p",
                0,
                1,
                errors
            );
            params.symmetricPhase = finiteNumber(
                rawParams.symmetricPhase,
                0,
                "symmetric-phase",
                "GHZ phase",
                -1000000,
                1000000,
                errors
            );
            totalDimension = 32;
            dimensions = [4, 8];
            stateKind = "mixed";
        } else if (family === "tiles_upb") {
            totalDimension = 9;
            dimensions = [3, 3];
            stateKind = "mixed";
        } else if (family === "random_density") {
            params.randomDimA = finiteInteger(
                rawParams.randomDimA,
                2,
                "random-dim-a",
                "Random-state subsystem A dimension",
                2,
                8,
                errors
            );
            params.randomDimB = finiteInteger(
                rawParams.randomDimB,
                3,
                "random-dim-b",
                "Random-state subsystem B dimension",
                2,
                8,
                errors
            );
            totalDimension = params.randomDimA * params.randomDimB;
            dimensions = [params.randomDimA, params.randomDimB];
            if (totalDimension > 64) {
                addError(
                    errors,
                    "random_dimension_limit",
                    "random-dim-b",
                    "The seeded random density matrix is capped at total dimension 64."
                );
            }
            params.randomRank = finiteInteger(
                rawParams.randomRank,
                totalDimension,
                "random-rank",
                "Random-state rank bound",
                1,
                totalDimension,
                errors
            );
            params.randomSeed = finiteInteger(
                rawParams.randomSeed,
                20260731,
                "random-seed",
                "Xoshiro seed",
                0,
                4294967295,
                errors
            );
            params.randomDistribution = enumValue(
                rawParams.randomDistribution,
                "hilbert_schmidt",
                ["hilbert_schmidt", "bures"],
                "random-distribution",
                "Random-state distribution",
                errors
            );
            params.randomReal = booleanSetting(
                rawParams,
                "randomReal",
                false,
                "random-real",
                "Real random-state output",
                errors
            );
            stateKind = "mixed";
        }

        if (totalDimension > 256 && family !== "ghz") {
            addError(
                errors,
                "resource_limit",
                "state-family",
                "The generated dense analysis is capped at total dimension 256."
            );
        }

        var rawAnalysis = objectOrEmpty(source.analysis);
        var analysis = {};
        var selectedCount = 0;
        for (var analysisIndex = 0; analysisIndex < ANALYSIS_ORDER.length; analysisIndex += 1) {
            var key = ANALYSIS_ORDER[analysisIndex];
            analysis[key] = booleanValue(rawAnalysis, key, key === "pipeline", errors);
            if (analysis[key]) {
                selectedCount += 1;
            }
        }
        analysis.ppt_witness = booleanValue(rawAnalysis, "ppt_witness", false, errors);
        if (analysis.ppt_witness && !analysis.ppt) {
            addError(
                errors,
                "witness_requires_ppt",
                "analysis-ppt-witness",
                "Select the PPT criterion before requesting its decomposable witness."
            );
        }
        if (selectedCount === 0) {
            addError(
                errors,
                "analysis_required",
                "analysis-pipeline",
                "Select at least one analysis or measure."
            );
        }
        if (analysis.schmidt && stateKind !== "pure") {
            addError(
                errors,
                "schmidt_requires_pure",
                "analysis-schmidt",
                "Schmidt decomposition is available only for a generated pure state."
            );
        }
        if (analysis.separability && totalDimension > 64) {
            addError(
                errors,
                "separability_resource_limit",
                "analysis-separability",
                "The composite separability search is capped at total dimension 64."
            );
        }

        var rawOptions = objectOrEmpty(source.options);
        var options = {
            separabilityProfile: enumValue(
                rawOptions.separabilityProfile,
                "balanced_core",
                ["ppt_only", "fast_detection", "balanced_core", "separable_ball"],
                "analysis-separability-profile",
                "Separability profile",
                errors
            ),
        };

        var rawOutput = objectOrEmpty(source.output);
        var output = {
            result_helpers: booleanSetting(
                rawOutput,
                "result_helpers",
                false,
                "output-result-helpers",
                "Result-helper output",
                errors
            ),
            include_assertions: booleanSetting(
                rawOutput,
                "include_assertions",
                false,
                "output-include-assertions",
                "Structural assertions",
                errors
            ),
        };

        var atol = finiteNumber(
            source.atol,
            1e-12,
            "analysis-atol",
            "Absolute tolerance",
            0,
            0.1,
            errors
        );
        var rtol = finiteNumber(
            source.rtol,
            0,
            "analysis-rtol",
            "Relative tolerance",
            0,
            0.1,
            errors
        );

        var config = {
            version: SCHEMA_VERSION,
            family: family,
            params: params,
            analysis: analysis,
            atol: atol,
            rtol: rtol,
            options: options,
            output: output,
            derived: {
                dimensions: dimensions,
                totalDimension: totalDimension,
                stateKind: stateKind,
            },
        };

        return { ok: errors.length === 0, errors: errors, config: config };
    }

    function portableConfig(config) {
        var analysis = {};
        for (var index = 0; index < ANALYSIS_ORDER.length; index += 1) {
            analysis[ANALYSIS_ORDER[index]] = config.analysis[ANALYSIS_ORDER[index]];
        }
        analysis.ppt_witness = config.analysis.ppt_witness;
        return {
            version: SCHEMA_VERSION,
            family: config.family,
            params: clone(config.params),
            analysis: analysis,
            atol: config.atol,
            rtol: config.rtol,
            options: {
                separabilityProfile: config.options.separabilityProfile,
            },
            output: {
                result_helpers: config.output.result_helpers,
                include_assertions: config.output.include_assertions,
            },
        };
    }

    function serializeConfig(raw) {
        var normalized = normalizeConfig(raw);
        if (!normalized.ok) {
            return {
                ok: false,
                errors: normalized.errors,
                config: normalized.config,
                text: "",
            };
        }
        var portable = portableConfig(normalized.config);
        var text = JSON.stringify(portable, null, 2) + "\n";
        if (utf8ByteLength(text) > MAX_CONFIG_BYTES) {
            return {
                ok: false,
                errors: [
                    {
                        code: "config_too_large",
                        field: "config-json",
                        message: "The canonical configuration exceeds the 16 KiB limit.",
                    },
                ],
                config: normalized.config,
                text: "",
            };
        }
        return { ok: true, errors: [], config: normalized.config, text: text };
    }

    function validateKnownKeys(object, allowed, path, errors) {
        if (!object || typeof object !== "object" || Array.isArray(object)) {
            addError(
                errors,
                "invalid_config_shape",
                "config-json",
                path + " must be a JSON object."
            );
            return;
        }
        for (var key in object) {
            if (hasOwn(object, key) && allowed.indexOf(key) === -1) {
                addError(
                    errors,
                    "unknown_config_key",
                    "config-json",
                    "Unsupported configuration key " + path + "." + key + "."
                );
            }
        }
    }

    function validateParsedConfig(source, errors) {
        validateKnownKeys(
            source,
            ["version", "family", "params", "analysis", "atol", "rtol", "options", "output"],
            "config",
            errors
        );
        if (hasOwn(source, "params")) {
            var allowedParams =
                typeof source.family === "string" &&
                hasOwn(FAMILY_PARAM_KEYS, source.family)
                    ? FAMILY_PARAM_KEYS[source.family]
                    : [].concat.apply(
                          [],
                          FAMILY_ORDER.map(function (family) {
                              return FAMILY_PARAM_KEYS[family];
                          })
                      );
            validateKnownKeys(
                source.params,
                allowedParams,
                "config.params",
                errors
            );
        }
        if (hasOwn(source, "analysis")) {
            validateKnownKeys(
                source.analysis,
                ANALYSIS_ORDER.concat(["ppt_witness"]),
                "config.analysis",
                errors
            );
        }
        if (hasOwn(source, "options")) {
            validateKnownKeys(
                source.options,
                ["separabilityProfile"],
                "config.options",
                errors
            );
        }
        if (hasOwn(source, "output")) {
            validateKnownKeys(
                source.output,
                ["result_helpers", "include_assertions"],
                "config.output",
                errors
            );
        }
    }

    function parseConfig(text) {
        var errors = [];
        if (typeof text !== "string") {
            addError(
                errors,
                "invalid_json",
                "config-json",
                "The configuration must be supplied as JSON text."
            );
            return { ok: false, errors: errors, config: null };
        }
        if (utf8ByteLength(text) > MAX_CONFIG_BYTES) {
            addError(
                errors,
                "config_too_large",
                "config-json",
                "The configuration exceeds the 16 KiB limit."
            );
            return { ok: false, errors: errors, config: null };
        }

        var parsed;
        try {
            parsed = JSON.parse(text);
        } catch (error) {
            addError(
                errors,
                "invalid_json",
                "config-json",
                "The configuration is not valid JSON."
            );
            return { ok: false, errors: errors, config: null };
        }
        if (!parsed || typeof parsed !== "object" || Array.isArray(parsed)) {
            addError(
                errors,
                "invalid_json",
                "config-json",
                "The configuration JSON must contain one object."
            );
            return { ok: false, errors: errors, config: null };
        }

        validateParsedConfig(parsed, errors);
        var normalized = normalizeConfig(parsed);
        errors = errors.concat(normalized.errors);
        return {
            ok: errors.length === 0,
            errors: errors,
            config: normalized.config,
        };
    }

    function juliaFloat(number) {
        if (number === 0) {
            return "0.0";
        }
        if (Math.floor(number) === number) {
            return String(number) + ".0";
        }
        return String(number).replace("e+", "e");
    }

    function emitHeader(lines, config) {
        var family = config.family;
        lines.push("# QuantumEntanglementTools.jl generated example");
        lines.push("#");
        lines.push("# SPDX-FileCopyrightText: 2026 John Martin");
        lines.push("# SPDX-License-Identifier: BSD-3-Clause");
        lines.push("# Full template terms: https://github.com/aenictusGitHub/QuantumEntanglementTools/blob/main/LICENSE");
        lines.push("#");
        lines.push("# This deterministic template was generated in your browser.");
        lines.push("# It does not record or upload the selected parameters.");
        lines.push("# Substantial template text is covered by the license above.");
        lines.push("# You remain responsible for rights in parameters, comments, and data you add.");
        lines.push("#");
        lines.push("# Run in an environment where QuantumEntanglementTools is available:");
        lines.push("#   julia --startup-file=no --project=. qet_" + family + "_example.jl");
        lines.push("");
        lines.push("using LinearAlgebra");
        if (family === "random_density") {
            lines.push("using Random");
        }
        lines.push("using QuantumEntanglementTools");
        lines.push("");
    }

    function emitState(lines, config) {
        var p = config.params;
        var family = config.family;
        var pipelineInput = "rho";

        lines.push("# --- State construction -------------------------------------------------------");
        if (family === "product_basis") {
            lines.push("# Physics labels are zero-based; Julia array positions are one-based.");
            lines.push("dim_a, dim_b = " + p.dimA + ", " + p.dimB);
            lines.push("basis_label_a, basis_label_b = " + p.indexA + ", " + p.indexB);
            lines.push("ket_a = zeros(ComplexF64, dim_a)");
            lines.push("ket_b = zeros(ComplexF64, dim_b)");
            lines.push("ket_a[basis_label_a + 1] = 1");
            lines.push("ket_b[basis_label_b + 1] = 1");
            lines.push("psi = tensor_product(ket_a, ket_b)");
            lines.push("dims = (dim_a, dim_b)");
            lines.push("rho = psi * psi'");
            pipelineInput = "psi";
        } else if (family === "bell") {
            lines.push("bell_index = " + p.bellIndex);
            lines.push(
                "psi = bell_state(bell_index; normalized=true, sparse_output=false, T=Float64)"
            );
            lines.push("dims = (2, 2)");
            lines.push("rho = psi * psi'");
            pipelineInput = "psi";
        } else if (family === "diagonal_mixture") {
            lines.push("# Every term below is an explicit product-basis projector.");
            lines.push("local_dimension = " + p.mixtureDimension);
            lines.push(
                "local_basis = [ComplexF64[position == label ? 1 : 0 for position in 1:local_dimension] for label in 1:local_dimension]"
            );
            lines.push(
                "product_basis = [tensor_product(left, right) for left in local_basis for right in local_basis]"
            );
            if (p.mixtureDimension === 2) {
                lines.push("weights = [3 / 8, 1 / 8, 1 / 8, 3 / 8]");
            } else {
                lines.push("weights = [fill(1 / 8, 7); 1 / 16; 1 / 16]");
            }
            lines.push(
                "rho = sum(weight * (state * state') for (weight, state) in zip(weights, product_basis))"
            );
            lines.push("dims = (local_dimension, local_dimension)");
        } else if (family === "ghz") {
            lines.push("local_dimension = " + p.ghzDimension);
            lines.push("parties = " + p.ghzParties);
            lines.push("left_parties = " + p.ghzCut);
            lines.push(
                "psi = ghz_state(local_dimension, parties; sparse_output=false, T=Float64)"
            );
            lines.push(
                "dims = (local_dimension^left_parties, local_dimension^(parties - left_parties))"
            );
            lines.push("rho = psi * psi'");
            pipelineInput = "psi";
        } else if (family === "dicke") {
            lines.push("parties = " + p.dickeParties);
            lines.push("excitations = " + p.dickeExcitations);
            lines.push("left_parties = " + p.dickeCut);
            lines.push(
                "psi = dicke_state(parties, excitations; normalized=true, sparse_output=false, T=Float64)"
            );
            lines.push("dims = (2^left_parties, 2^(parties - left_parties))");
            lines.push("rho = psi * psi'");
            pipelineInput = "psi";
        } else if (family === "isotropic") {
            lines.push("local_dimension = " + p.isotropicDimension);
            lines.push("alpha = " + juliaFloat(p.isotropicAlpha));
            lines.push(
                "rho = isotropic_state(local_dimension, alpha; sparse_output=false)"
            );
            lines.push("dims = (local_dimension, local_dimension)");
        } else if (family === "werner") {
            lines.push("local_dimension = " + p.wernerDimension);
            lines.push("alpha = " + juliaFloat(p.wernerAlpha));
            lines.push("rho = werner_state(local_dimension, alpha; sparse_output=false)");
            lines.push("dims = (local_dimension, local_dimension)");
        } else if (family === "horodecki") {
            lines.push("a = " + juliaFloat(p.horodeckiA));
            lines.push(
                "dims = " +
                    (p.horodeckiDims === "3x3" ? "(3, 3)" : "(2, 4)")
            );
            lines.push("rho = horodecki_state(a; dims=dims)");
        } else if (family === "symmetric_sappt_ghz5") {
            lines.push("# Source: J. Louvet et al., \"Nonequivalence between absolute");
            lines.push("# separability and positive partial transposition in the symmetric");
            lines.push("# subspace,\" Phys. Rev. A 111, 042418 (2025),");
            lines.push("# https://doi.org/10.1103/PhysRevA.111.042418.");
            lines.push("# The state family and rounded W5 coefficients below are cited from");
            lines.push("# that paper; this template does not copy or rerun the source SDP.");
            lines.push("qubits = 5");
            lines.push("p = " + juliaFloat(p.symmetricP));
            lines.push("phase = " + juliaFloat(p.symmetricPhase));
            lines.push("ghz_dicke = zeros(ComplexF64, qubits + 1)");
            lines.push("ghz_dicke[1] = inv(sqrt(2))");
            lines.push("ghz_dicke[end] = cis(phase) / sqrt(2)");
            lines.push(
                "rho_dicke = p * Matrix{ComplexF64}(I, qubits + 1, qubits + 1) / (qubits + 1) +"
            );
            lines.push("            (1 - p) * (ghz_dicke * ghz_dicke')");
            lines.push(
                "dicke_basis = symmetric_subspace_basis(2, qubits; sparse_output=false)"
            );
            lines.push("rho = dicke_basis * rho_dicke * dicke_basis'");
            lines.push("# Group the first two qubits against the remaining three.");
            lines.push("dims = (4, 8)");
            lines.push("");
            lines.push("# Rounded, phase-matched symmetric-subspace witness from the paper.");
            lines.push("a_w, b_w, c_w = 0.0366656, -0.134595, -9.31947");
            lines.push(
                "w5 = Matrix(Diagonal(ComplexF64[a_w, b_w, 1, 1, b_w, a_w]))"
            );
            lines.push("w5[1, end] = c_w * cis(-phase)");
            lines.push("w5[end, 1] = conj(w5[1, end])");
            lines.push("w5_expectation = real(tr(w5 * rho_dicke))");
            lines.push('println("Published W5 expectation: ", w5_expectation)');
            lines.push(
                'println("  A negative value certifies entanglement; a nonnegative value is inconclusive.")'
            );
        } else if (family === "tiles_upb") {
            lines.push("# Exact integer local factors for the two-qutrit Tiles UPB.");
            lines.push("tiles_left = BigInt[");
            lines.push("    1 1 0 0 1");
            lines.push("    0 -1 0 1 1");
            lines.push("    0 0 1 -1 1");
            lines.push("]");
            lines.push("tiles_right = BigInt[");
            lines.push("    1 0 0 1 1");
            lines.push("    -1 0 1 0 1");
            lines.push("    0 1 -1 0 1");
            lines.push("]");
            lines.push(
                "tiles_upb_analysis = is_upb(tiles_left, tiles_right; normalization=:allow)"
            );
            lines.push("tiles_projector = zeros(Rational{BigInt}, 9, 9)");
            lines.push("for column in axes(tiles_left, 2)");
            lines.push(
                "    product_vector = Rational{BigInt}.(tensor_product(tiles_left[:, column], tiles_right[:, column]))"
            );
            lines.push(
                "    tiles_projector .+= (product_vector * product_vector') / dot(product_vector, product_vector)"
            );
            lines.push("end");
            lines.push(
                "tiles_complement = Matrix{Rational{BigInt}}(I, 9, 9) - tiles_projector"
            );
            lines.push("rho_exact = tiles_complement / tr(tiles_complement)");
            lines.push(
                "tiles_partial_transpose_exact = partial_transpose(rho_exact, (3, 3); systems=(2,))"
            );
            lines.push("tiles_exact_ppt = tiles_partial_transpose_exact == rho_exact");
            lines.push("rho = Float64.(rho_exact)");
            lines.push("dims = (3, 3)");
            lines.push(
                'println("Exact Tiles construction: trace=", tr(rho_exact), ", PPT=", tiles_exact_ppt, ", UPB status=", tiles_upb_analysis.status)'
            );
        } else if (family === "random_density") {
            lines.push("dims = (" + p.randomDimA + ", " + p.randomDimB + ")");
            lines.push("random_rank = " + p.randomRank);
            lines.push("random_seed = UInt64(" + p.randomSeed + ")");
            lines.push("random_distribution = :" + p.randomDistribution);
            lines.push("random_real = " + (p.randomReal ? "true" : "false"));
            lines.push(
                "rho = random_density_matrix(Xoshiro(random_seed), prod(dims); rank=random_rank, distribution=random_distribution, real=random_real)"
            );
            lines.push(
                'println("Seeded random density matrix: seed=", random_seed, ", distribution=", random_distribution, ", requested rank=", random_rank)'
            );
        }
        lines.push("ppt_systems = (2,)");
        lines.push("realignment_systems = (1,)");
        lines.push("");
        return pipelineInput;
    }

    function emitStructuralAssertions(lines, config) {
        if (!config.output.include_assertions) {
            return;
        }
        lines.push("# --- Structural assertions ----------------------------------------------------");
        lines.push("@assert length(dims) == 2");
        lines.push("@assert all(dimension -> dimension >= 2, dims)");
        lines.push("@assert size(rho) == (prod(dims), prod(dims))");
        lines.push("@assert all(isfinite, rho)");
        lines.push("@assert ishermitian(rho)");
        if (config.derived.stateKind === "pure") {
            lines.push("@assert length(psi) == prod(dims)");
            lines.push("@assert size(rho) == (length(psi), length(psi))");
        }
        if (config.family === "diagonal_mixture") {
            lines.push("@assert sum(weights) == 1");
        }
        if (config.family === "tiles_upb") {
            lines.push("@assert size(rho_exact) == (9, 9)");
            lines.push("@assert tr(rho_exact) == 1");
        }
        lines.push("");
    }

    function juliaSymbolTuple(values) {
        var symbols = [];
        for (var index = 0; index < values.length; index += 1) {
            symbols.push(":" + values[index]);
        }
        if (symbols.length === 1) {
            return "(" + symbols[0] + ",)";
        }
        return "(" + symbols.join(", ") + ")";
    }

    function emitResultHelper(lines, config) {
        if (!config.output.result_helpers) {
            return;
        }
        lines.push("# Keep the native status visible; the helpers interpret but never replace it.");
        lines.push("function print_result_summary(label, result)");
        lines.push('    println(label, ": raw status=", result.status)');
        lines.push(
            '    println("  conclusion=", conclusion(result), ", conclusive=", is_conclusive(result), ", certified=", is_certified(result))'
        );
        lines.push('    println("  ", explain(result))');
        lines.push("end");
        lines.push("");
    }

    function emitAnalysis(lines, config, pipelineInput) {
        var a = config.analysis;
        var atol = juliaFloat(config.atol);
        var rtol = juliaFloat(config.rtol);
        var toleranceKeywords = "; atol=" + atol + ", rtol=" + rtol;
        var includeAssertions = config.output.include_assertions;
        var helperOutput = config.output.result_helpers;

        lines.push("# --- Requested analysis -------------------------------------------------------");
        emitResultHelper(lines, config);

        if (a.validation) {
            lines.push(
                "validation_report = validate_density_matrix(rho, dims; atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ", allow_densify=false, max_dense_entries=65_536)"
            );
            lines.push(
                'println("Validation: raw status=", validation_report.status, ", valid=", validation_report.valid, ", complete=", validation_report.complete)'
            );
            lines.push(
                'println("  trace=", validation_report.trace_value, ", minimum eigenvalue=", validation_report.minimum_eigenvalue, ", spectral analysis=", validation_report.spectral_analysis)'
            );
            lines.push("for message in validation_report.messages");
            lines.push('    println("  - ", message)');
            lines.push("end");
            if (includeAssertions) {
                lines.push("@assert validation_report.valid");
            }
            lines.push("");
        }

        if (a.marginals) {
            lines.push("marginal_a = partial_trace(rho, dims; trace_out=(2,))");
            lines.push("marginal_b = partial_trace(rho, dims; trace_out=(1,))");
            lines.push(
                'println("Marginal A: size=", size(marginal_a), ", trace=", tr(marginal_a))'
            );
            lines.push(
                'println("Marginal B: size=", size(marginal_b), ", trace=", tr(marginal_b))'
            );
            if (includeAssertions) {
                lines.push("@assert size(marginal_a) == (dims[1], dims[1])");
                lines.push("@assert size(marginal_b) == (dims[2], dims[2])");
                lines.push("@assert all(isfinite, marginal_a)");
                lines.push("@assert all(isfinite, marginal_b)");
            }
            lines.push("");
        }

        if (a.schmidt) {
            lines.push(
                "schmidt_result = schmidt_decomposition(psi, dims; allow_densify=false)"
            );
            lines.push(
                "schmidt_numeric_rank = schmidt_rank(psi, dims; atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ", allow_densify=false)"
            );
            lines.push(
                'println("Schmidt coefficients: ", schmidt_result.coefficients)'
            );
            lines.push('println("Schmidt numerical rank: ", schmidt_numeric_rank)');
            if (includeAssertions) {
                lines.push(
                    "@assert length(schmidt_result.coefficients) == min(dims...)"
                );
            }
            lines.push("");
        }

        if (a.backend_status) {
            lines.push("backend_readiness = backend_status()");
            lines.push(
                'println("Native backend: loaded=", backend_readiness.entanglement.native.loaded, ", ready=", backend_readiness.entanglement.native.ready)'
            );
            lines.push(
                'println("  ", backend_readiness.entanglement.native.message)'
            );
            lines.push(
                'println("EntanglementDetection backend: loaded=", backend_readiness.entanglement.entanglement_detection.loaded, ", ready=", backend_readiness.entanglement.entanglement_detection.ready)'
            );
            lines.push(
                'println("  ", backend_readiness.entanglement.entanglement_detection.message)'
            );
            lines.push(
                'println("Optimization backend: configured=", backend_readiness.optimization.configured, ", ready=", backend_readiness.optimization.ready)'
            );
            lines.push('println("  ", backend_readiness.optimization.message)');
            if (includeAssertions) {
                lines.push("@assert backend_readiness.entanglement.native.loaded");
                lines.push("@assert backend_readiness.entanglement.native.ready");
                lines.push("@assert !backend_readiness.optimization.configured");
            }
            lines.push("");
        }

        if (a.pipeline) {
            if (config.derived.stateKind === "pure") {
                lines.push(
                    "report = analyze_entanglement(" +
                        pipelineInput +
                        ", dims" +
                        toleranceKeywords +
                        ")"
                );
            } else {
                lines.push(
                    "report = analyze_entanglement(rho, dims; systems=ppt_systems, atol=" +
                        atol +
                        ", rtol=" +
                        rtol +
                        ")"
                );
            }
            if (helperOutput) {
                lines.push('print_result_summary("Pipeline", report)');
                lines.push(
                    'println("  certificate kind=", report.certificate_kind)'
                );
            } else {
                lines.push(
                    'println("Pipeline: status=", report.status, ", certified=", report.certified, ", certificate=", report.certificate_kind)'
                );
                lines.push('println("  ", report.message)');
            }
            lines.push("for attempt in report.attempts");
            lines.push(
                '    println("  - ", attempt.method, ": ", attempt.status, " (certified=", attempt.certified, ")")'
            );
            lines.push("end");
            if (includeAssertions) {
                lines.push("@assert report isa EntanglementReport");
            }
            lines.push("");
        }

        if (a.ppt) {
            lines.push(
                "ppt_result = ppt_criterion(rho, dims; systems=ppt_systems, atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ")"
            );
            if (helperOutput) {
                lines.push('print_result_summary("PPT", ppt_result)');
                lines.push(
                    'println("  minimum eigenvalue=", ppt_result.value, ", tolerance=", ppt_result.tolerance)'
                );
            } else {
                lines.push(
                    'println("PPT: status=", ppt_result.status, ", minimum eigenvalue=", ppt_result.value, ", tolerance=", ppt_result.tolerance)'
                );
                lines.push('println("  ", ppt_result.message)');
            }
            lines.push("if ppt_result.status === CriterionSatisfied");
            lines.push(
                '    println("  Passing PPT alone is not a separability certificate in general.")'
            );
            lines.push("end");
            if (includeAssertions) {
                lines.push("@assert ppt_result isa CriterionResult");
            }
            lines.push("");
        }

        if (a.realignment) {
            lines.push(
                "realignment_result = realignment_criterion(rho, dims; systems=realignment_systems, atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ")"
            );
            if (helperOutput) {
                lines.push(
                    'print_result_summary("Realignment", realignment_result)'
                );
                lines.push(
                    'println("  trace norm=", realignment_result.value, ", boundary=", realignment_result.threshold)'
                );
            } else {
                lines.push(
                    'println("Realignment: status=", realignment_result.status, ", trace norm=", realignment_result.value, ", boundary=", realignment_result.threshold)'
                );
                lines.push('println("  ", realignment_result.message)');
            }
            if (includeAssertions) {
                lines.push("@assert realignment_result isa CriterionResult");
            }
            lines.push("");
        }

        if (a.reduction) {
            lines.push(
                "reduction_result = reduction_criterion(rho, dims; side=:both, atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ")"
            );
            if (helperOutput) {
                lines.push('print_result_summary("Reduction", reduction_result)');
                lines.push(
                    'println("  minimum eigenvalue=", reduction_result.value)'
                );
            } else {
                lines.push(
                    'println("Reduction: status=", reduction_result.status, ", minimum eigenvalue=", reduction_result.value)'
                );
                lines.push('println("  ", reduction_result.message)');
            }
            if (includeAssertions) {
                lines.push("@assert reduction_result isa CriterionResult");
            }
            lines.push("");
        }

        if (a.separable_ball) {
            lines.push(
                "ball_result = in_separable_ball(rho, dims; atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ")"
            );
            if (helperOutput) {
                lines.push('print_result_summary("Separable ball", ball_result)');
                lines.push(
                    'println("  purity=", ball_result.purity, ", boundary=", ball_result.boundary)'
                );
            } else {
                lines.push(
                    'println("Separable ball: status=", ball_result.status, ", purity=", ball_result.purity, ", boundary=", ball_result.boundary)'
                );
                lines.push('println("  ", ball_result.message)');
            }
            lines.push(
                'println("  Only :separable_certified is conclusive; :outside_ball is not entanglement.")'
            );
            if (includeAssertions) {
                lines.push("@assert ball_result isa SeparableBallResult");
            }
            lines.push("");
        }

        if (a.separability) {
            var profileStrategies =
                SEPARABILITY_PROFILES[config.options.separabilityProfile];
            var strategyTuple = juliaSymbolTuple(profileStrategies);
            lines.push(
                "# Explicit dependency-free strategy profile: " +
                    config.options.separabilityProfile
            );
            lines.push("separability_strategies = " + strategyTuple);
            lines.push(
                "available_strategies = available_separability_strategies()"
            );
            lines.push("for strategy in separability_strategies");
            lines.push(
                '    strategy in available_strategies || error("Generated strategy is no longer available: $(strategy)")'
            );
            lines.push("    strategy_metadata = describe_strategy(strategy)");
            lines.push(
                '    strategy_metadata.optional_dependency && error("Generated core profile unexpectedly requires an optional dependency: $(strategy)")'
            );
            lines.push(
                '    strategy_metadata.rng_required && error("Generated deterministic profile unexpectedly requires an RNG: $(strategy)")'
            );
            lines.push(
                '    println("Strategy ", strategy, ": cost=", strategy_metadata.cost, ", certificate directions=", strategy_metadata.certificate_directions)'
            );
            lines.push("end");
            lines.push(
                "separability_report = if iszero(imag(tr(rho))) && real(tr(rho)) == 1"
            );
            lines.push("    try");
            lines.push(
                "        is_separable(rho, dims; strategies=separability_strategies, atol=" +
                    atol +
                    ", rtol=" +
                    rtol +
                    ", allow_densify=false, max_dense_entries=65_536, max_work=10_000_000)"
            );
            lines.push("    catch error");
            lines.push("        error isa InterruptException && rethrow()");
            lines.push(
                '        println("Separability search unavailable within the selected conservative limits; the state was not normalized or repaired: ", sprint(showerror, error))'
            );
            lines.push("        nothing");
            lines.push("    end");
            lines.push("else");
            lines.push(
                '    println("Separability search skipped: is_separable requires an exactly represented real unit trace; the input was not normalized or repaired.")'
            );
            lines.push("    nothing");
            lines.push("end");
            lines.push("if separability_report !== nothing");
            if (helperOutput) {
                lines.push(
                    '    print_result_summary("Separability search", separability_report)'
                );
                lines.push(
                    '    println("  certificate kind=", separability_report.certificate_kind)'
                );
            } else {
                lines.push(
                    '    println("Separability search: status=", separability_report.status, ", certified=", separability_report.certified, ", certificate=", separability_report.certificate_kind)'
                );
                lines.push('    println("  ", separability_report.message)');
            }
            lines.push("    for attempt in separability_report.attempts");
            lines.push(
                '        println("  - ", attempt.method, ": ", attempt.status, " (certified=", attempt.certified, ")")'
            );
            lines.push("    end");
            if (includeAssertions) {
                lines.push(
                    "    @assert separability_report isa EntanglementReport"
                );
            }
            lines.push("end");
            lines.push("");
        }

        if (a.measures) {
            lines.push("function evaluate_without_repair(label, operation)");
            lines.push("    try");
            lines.push("        return operation()");
            lines.push("    catch error");
            lines.push("        error isa InterruptException && rethrow()");
            lines.push(
                '        println(label, " unavailable without altering the state: ", sprint(showerror, error))'
            );
            lines.push("        return nothing");
            lines.push("    end");
            lines.push("end");
            lines.push(
                'state_purity = evaluate_without_repair("Purity", () -> purity(rho; atol=' +
                    atol +
                    ", rtol=" +
                    rtol +
                    "))"
            );
            lines.push(
                'entropy_bits = evaluate_without_repair("Entropy", () -> von_neumann_entropy(rho; base=2, atol=' +
                    atol +
                    ", rtol=" +
                    rtol +
                    "))"
            );
            lines.push(
                'state_negativity = evaluate_without_repair("Negativity", () -> negativity(rho, dims; systems=ppt_systems, atol=' +
                    atol +
                    ", rtol=" +
                    rtol +
                    "))"
            );
            lines.push(
                'log_negativity = evaluate_without_repair("Logarithmic negativity", () -> logarithmic_negativity(rho, dims; systems=ppt_systems, base=2, atol=' +
                    atol +
                    ", rtol=" +
                    rtol +
                    "))"
            );
            lines.push(
                'println("Measures: purity=", state_purity, ", entropy(base 2)=", entropy_bits)'
            );
            lines.push(
                'println("  negativity=", state_negativity, ", logarithmic negativity(base 2)=", log_negativity)'
            );
            lines.push(
                'println("  A zero negativity is not, by itself, a separability certificate.")'
            );
            lines.push("");
        }

        if (a.ppt_witness) {
            lines.push("# Construct the decomposable witness W = (|eta><eta|)^(T_systems).");
            lines.push("ppt_witness = nothing");
            lines.push("if ppt_result.status === CriterionEntanglementDetected");
            lines.push("    eta = ppt_result.witness");
            lines.push("    positive_seed = eta * eta'");
            lines.push(
                "    ppt_witness = partial_transpose(positive_seed, dims; systems=ppt_systems)"
            );
            lines.push("    witness_expectation = real(tr(ppt_witness * rho))");
            lines.push(
                '    println("Decomposable PPT-witness expectation: ", witness_expectation)'
            );
            if (includeAssertions) {
                lines.push("    @assert size(ppt_witness) == size(rho)");
            }
            lines.push("else");
            lines.push(
                '    println("No decomposable PPT witness: PPT did not robustly detect entanglement.")'
            );
            lines.push("end");
            lines.push("");
        }
    }

    function noticesFor(config) {
        var notices = [];
        if (config.analysis.ppt || config.analysis.realignment || config.analysis.reduction) {
            notices.push(
                "PPT, realignment, and reduction are one-sided tests: a robust violation certifies entanglement, but passing does not generally certify separability."
            );
        }
        if (config.analysis.separable_ball) {
            notices.push(
                "The separable-ball result is conclusive only when its status is :separable_certified; :outside_ball and :unknown remain inconclusive."
            );
        }
        if (config.analysis.measures) {
            notices.push(
                "Negativity and logarithmic negativity are quantitative diagnostics, not complete separability tests."
            );
        }
        if (config.analysis.validation) {
            notices.push(
                "Validation reports the supplied matrix without normalizing, symmetrizing, clipping, or otherwise repairing it."
            );
        }
        if (config.analysis.separability) {
            notices.push(
                "The composite separability search uses an explicit dependency-free strategy tuple and fixed resource limits. An unavailable or inconclusive route remains unknown."
            );
        }
        if (config.analysis.schmidt) {
            notices.push(
                "Schmidt diagnostics are generated only for pure-state families and use the selected bipartition."
            );
        }
        if (config.analysis.backend_status) {
            notices.push(
                "Backend readiness describes software availability only; it is not a mathematical conclusion about the state."
            );
        }
        if (config.family === "tiles_upb") {
            notices.push(
                "The Tiles complement is constructed first in exact BigInt/Rational arithmetic and only then copied to Float64 for numerical analyses."
            );
        }
        if (config.family === "random_density") {
            notices.push(
                "The random density matrix uses a local Xoshiro instance with the displayed seed and never mutates Julia's default random stream."
            );
        }
        if (config.output.result_helpers) {
            notices.push(
                "Result helpers display each native status before its conservative conclusion, conclusiveness, certification, and explanation."
            );
        }
        if (config.family === "symmetric_sappt_ghz5") {
            notices.push(
                "The W5 coefficients are the rounded published symmetric-subspace witness. The generated code evaluates it but does not rerun the source SDP or prove SAPPT numerically."
            );
        }
        if (config.derived.stateKind === "pure" && config.analysis.pipeline) {
            notices.push(
                "The pipeline receives the pure vector and uses its Schmidt structure; matrix-only criteria receive rho = psi * psi'."
            );
        }
        return notices;
    }

    function generate(raw) {
        var normalized = normalizeConfig(raw);
        if (!normalized.ok) {
            return {
                ok: false,
                errors: normalized.errors,
                config: normalized.config,
                code: "",
                filename: "",
                summary: null,
                notices: [],
            };
        }

        var config = normalized.config;
        var lines = [];
        emitHeader(lines, config);
        var pipelineInput = emitState(lines, config);
        emitStructuralAssertions(lines, config);
        emitAnalysis(lines, config, pipelineInput);
        lines.push("# End of generated example.");

        var selectedAnalyses = [];
        for (var index = 0; index < ANALYSIS_ORDER.length; index += 1) {
            var key = ANALYSIS_ORDER[index];
            if (config.analysis[key]) {
                selectedAnalyses.push(ANALYSIS_LABELS[key]);
            }
        }
        if (config.analysis.ppt_witness) {
            selectedAnalyses.push("decomposable PPT witness");
        }

        return {
            ok: true,
            errors: [],
            config: config,
            code: lines.join("\n") + "\n",
            filename: "qet_" + config.family + "_example.jl",
            summary: {
                family: FAMILY_LABELS[config.family],
                stateKind: config.derived.stateKind,
                dimensions: config.derived.dimensions.slice(),
                totalDimension: config.derived.totalDimension,
                analyses: selectedAnalyses,
            },
            notices: noticesFor(config),
        };
    }

    function presetConfig(id) {
        if (!hasOwn(PRESETS, id)) {
            return null;
        }
        return clone(PRESETS[id].config);
    }

    function presetList() {
        var result = [];
        for (var id in PRESETS) {
            if (hasOwn(PRESETS, id)) {
                result.push({
                    id: id,
                    label: PRESETS[id].label,
                    description: PRESETS[id].description,
                });
            }
        }
        return result;
    }

    function familyList() {
        return FAMILY_ORDER.map(function (id) {
            return { id: id, label: FAMILY_LABELS[id] };
        });
    }

    return {
        schemaVersion: SCHEMA_VERSION,
        normalizeConfig: normalizeConfig,
        serializeConfig: serializeConfig,
        parseConfig: parseConfig,
        generate: generate,
        presetConfig: presetConfig,
        presetList: presetList,
        familyList: familyList,
    };
});

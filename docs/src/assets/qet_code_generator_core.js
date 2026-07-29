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

    var SCHEMA_VERSION = 1;
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
    };
    var ANALYSIS_ORDER = [
        "pipeline",
        "ppt",
        "realignment",
        "reduction",
        "separable_ball",
        "measures",
    ];
    var ANALYSIS_LABELS = {
        pipeline: "certificate-aware pipeline",
        ppt: "PPT criterion",
        realignment: "realignment / CCNR",
        reduction: "reduction criterion",
        separable_ball: "separable-ball certificate",
        measures: "entropy and negativities",
    };

    function clone(value) {
        return JSON.parse(JSON.stringify(value));
    }

    function hasOwn(object, key) {
        return Object.prototype.hasOwnProperty.call(object, key);
    }

    function objectOrEmpty(value) {
        return value && typeof value === "object" && !Array.isArray(value) ? value : {};
    }

    function rawConfig(family, params, analysis, atol, rtol) {
        return {
            version: SCHEMA_VERSION,
            family: family,
            params: params || {},
            analysis: analysis || {},
            atol: atol === undefined ? 1e-12 : atol,
            rtol: rtol === undefined ? 0 : rtol,
        };
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
                    ppt_witness: false,
                }
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
                    ppt_witness: true,
                }
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
                    ppt_witness: false,
                }
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
                    ppt_witness: false,
                },
                1e-12,
                1e-10
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
                    ppt_witness: false,
                },
                1e-12,
                0
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

    function booleanValue(source, key, fallback, errors) {
        if (!hasOwn(source, key)) {
            return fallback;
        }
        if (typeof source[key] !== "boolean") {
            addError(
                errors,
                "boolean_required",
                "analysis-" + key.replace(/_/g, "-"),
                (ANALYSIS_LABELS[key] || key) +
                    " must be selected or cleared explicitly."
            );
            return fallback;
        }
        return source[key];
    }

    function normalizeConfig(raw) {
        var errors = [];
        var source = objectOrEmpty(raw);
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

        if (
            hasOwn(source, "version") &&
            source.version !== undefined &&
            source.version !== SCHEMA_VERSION
        ) {
            addError(
                errors,
                "unsupported_schema",
                "state-family",
                "This configuration uses an unsupported generator schema version."
            );
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
            derived: {
                dimensions: dimensions,
                totalDimension: totalDimension,
                stateKind: stateKind,
            },
        };

        return { ok: errors.length === 0, errors: errors, config: config };
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

    function emitHeader(lines, family) {
        lines.push("# QuantumEntanglementTools.jl generated example");
        lines.push("#");
        lines.push("# This deterministic template was generated in your browser.");
        lines.push("# It does not record or upload the selected parameters.");
        lines.push("# Substantial template text is covered by the repository's BSD 3-Clause license.");
        lines.push("#");
        lines.push("# Run in an environment where QuantumEntanglementTools is available:");
        lines.push("#   julia --startup-file=no --project=. qet_" + family + "_example.jl");
        lines.push("");
        lines.push("using LinearAlgebra");
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
            lines.push("@assert sum(weights) == 1");
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
            lines.push("# Phys. Rev. A 111, 042418 (2025), five-qubit GHZ representative.");
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
        }
        lines.push("ppt_systems = (2,)");
        lines.push("realignment_systems = (1,)");
        lines.push("");
        return pipelineInput;
    }

    function emitAnalysis(lines, config, pipelineInput) {
        var a = config.analysis;
        var toleranceKeywords =
            "; atol=" + juliaFloat(config.atol) + ", rtol=" + juliaFloat(config.rtol);

        lines.push("# --- Requested analysis -------------------------------------------------------");
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
                        juliaFloat(config.atol) +
                        ", rtol=" +
                        juliaFloat(config.rtol) +
                        ")"
                );
            }
            lines.push(
                'println("Pipeline: status=", report.status, ", certified=", report.certified, ", certificate=", report.certificate_kind)'
            );
            lines.push('println("  ", report.message)');
            lines.push("for attempt in report.attempts");
            lines.push(
                '    println("  - ", attempt.method, ": ", attempt.status, " (certified=", attempt.certified, ")")'
            );
            lines.push("end");
            lines.push("");
        }

        if (a.ppt) {
            lines.push(
                "ppt_result = ppt_criterion(rho, dims; systems=ppt_systems, atol=" +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    ")"
            );
            lines.push(
                'println("PPT: status=", ppt_result.status, ", minimum eigenvalue=", ppt_result.value, ", tolerance=", ppt_result.tolerance)'
            );
            lines.push('println("  ", ppt_result.message)');
            lines.push("if ppt_result.status === CriterionSatisfied");
            lines.push(
                '    println("  Passing PPT alone is not a separability certificate in general.")'
            );
            lines.push("end");
            lines.push("");
        }

        if (a.realignment) {
            lines.push(
                "realignment_result = realignment_criterion(rho, dims; systems=realignment_systems, atol=" +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    ")"
            );
            lines.push(
                'println("Realignment: status=", realignment_result.status, ", trace norm=", realignment_result.value, ", boundary=", realignment_result.threshold)'
            );
            lines.push('println("  ", realignment_result.message)');
            lines.push("");
        }

        if (a.reduction) {
            lines.push(
                "reduction_result = reduction_criterion(rho, dims; side=:both, atol=" +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    ")"
            );
            lines.push(
                'println("Reduction: status=", reduction_result.status, ", minimum eigenvalue=", reduction_result.value)'
            );
            lines.push('println("  ", reduction_result.message)');
            lines.push("");
        }

        if (a.separable_ball) {
            lines.push(
                "ball_result = in_separable_ball(rho, dims; atol=" +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    ")"
            );
            lines.push(
                'println("Separable ball: status=", ball_result.status, ", purity=", ball_result.purity, ", boundary=", ball_result.boundary)'
            );
            lines.push('println("  ", ball_result.message)');
            lines.push(
                'println("  Only :separable_certified is conclusive; :outside_ball is not entanglement.")'
            );
            lines.push("");
        }

        if (a.measures) {
            lines.push("function evaluate_without_repair(label, operation)");
            lines.push("    try");
            lines.push("        return operation()");
            lines.push("    catch error");
            lines.push(
                '        println(label, " unavailable without altering the state: ", sprint(showerror, error))'
            );
            lines.push("        return nothing");
            lines.push("    end");
            lines.push("end");
            lines.push(
                'state_purity = evaluate_without_repair("Purity", () -> purity(rho; atol=' +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    "))"
            );
            lines.push(
                'entropy_bits = evaluate_without_repair("Entropy", () -> von_neumann_entropy(rho; base=2, atol=' +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    "))"
            );
            lines.push(
                'state_negativity = evaluate_without_repair("Negativity", () -> negativity(rho, dims; systems=ppt_systems, atol=' +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
                    "))"
            );
            lines.push(
                'log_negativity = evaluate_without_repair("Logarithmic negativity", () -> logarithmic_negativity(rho, dims; systems=ppt_systems, base=2, atol=' +
                    juliaFloat(config.atol) +
                    ", rtol=" +
                    juliaFloat(config.rtol) +
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
            lines.push("    @assert witness_expectation < -ppt_result.tolerance");
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
        emitHeader(lines, config.family);
        var pipelineInput = emitState(lines, config);
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
        generate: generate,
        presetConfig: presetConfig,
        presetList: presetList,
        familyList: familyList,
    };
});

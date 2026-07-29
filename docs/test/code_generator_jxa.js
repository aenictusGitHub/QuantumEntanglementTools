ObjC.import("Foundation");

function readUtf8(path) {
    var error = Ref();
    var value = $.NSString.stringWithContentsOfFileEncodingError(
        path,
        $.NSUTF8StringEncoding,
        error
    );
    if (!value) {
        throw new Error("Could not read " + path + ": " + ObjC.unwrap(error[0]));
    }
    return ObjC.unwrap(value);
}

function writeUtf8(path, value) {
    var error = Ref();
    var written = $(value).writeToFileAtomicallyEncodingError(
        path,
        true,
        $.NSUTF8StringEncoding,
        error
    );
    if (!written) {
        throw new Error("Could not write " + path + ": " + ObjC.unwrap(error[0]));
    }
}

function run(arguments) {
    var repositoryRoot = ObjC.unwrap($.NSFileManager.defaultManager.currentDirectoryPath);
    var outputDirectory =
        arguments.length > 0 ? String(arguments[0]) : "/tmp/qet-generator-smoke";
    var corePath =
        repositoryRoot + "/docs/src/assets/qet_code_generator_core.js";
    var uiPath =
        repositoryRoot + "/docs/src/assets/qet_code_generator_ui.js";
    var casesPath = repositoryRoot + "/docs/test/code_generator_cases.js";
    var pagePath = repositoryRoot + "/docs/src/code_generator.md";

    var globalEval = eval;
    globalEval(readUtf8(corePath));
    globalEval(readUtf8(casesPath));
    new Function(readUtf8(uiPath));

    var result = QETCodeGeneratorCases.run(QETCodeGenerator);
    var pageSource = readUtf8(pagePath);
    var uiSource = readUtf8(uiPath);
    var documentationFailures = QETCodeGeneratorCases.validateDocumentation(
        pageSource,
        uiSource
    );
    for (
        var failureIndex = 0;
        failureIndex < documentationFailures.length;
        failureIndex += 1
    ) {
        result.failures.push(documentationFailures[failureIndex]);
    }
    if (result.failures.length > 0) {
        throw new Error("Code-generator checks failed:\n" + result.failures.join("\n"));
    }

    var manager = $.NSFileManager.defaultManager;
    var directoryError = Ref();
    var created = manager.createDirectoryAtPathWithIntermediateDirectoriesAttributesError(
        outputDirectory,
        true,
        undefined,
        directoryError
    );
    if (!created) {
        throw new Error(
            "Could not create " +
                outputDirectory +
                ": " +
                ObjC.unwrap(directoryError[0])
        );
    }

    var outputPath = outputDirectory + "/generated_smoke.jl";
    writeUtf8(outputPath, result.bundle);
    return (
        "Code-generator core checks passed: " +
        result.checks +
        " assertions. Wrote " +
        outputPath
    );
}

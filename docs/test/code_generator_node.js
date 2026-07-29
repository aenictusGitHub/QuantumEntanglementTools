"use strict";

var fs = require("fs");
var path = require("path");

var repositoryRoot = path.resolve(__dirname, "..", "..");
var corePath = path.join(
    repositoryRoot,
    "docs",
    "src",
    "assets",
    "qet_code_generator_core.js"
);
var uiPath = path.join(
    repositoryRoot,
    "docs",
    "src",
    "assets",
    "qet_code_generator_ui.js"
);
var pagePath = path.join(repositoryRoot, "docs", "src", "code_generator.md");
var casesPath = path.join(__dirname, "code_generator_cases.js");

var core = require(corePath);
var cases = require(casesPath);
var uiSource = fs.readFileSync(uiPath, "utf8");
var pageSource = fs.readFileSync(pagePath, "utf8");

// Compile the thin DOM adapter without executing it.
new Function(uiSource);

var result = cases.run(core);
var documentationFailures = cases.validateDocumentation(pageSource, uiSource);
for (
    var failureIndex = 0;
    failureIndex < documentationFailures.length;
    failureIndex += 1
) {
    result.failures.push(documentationFailures[failureIndex]);
}

if (result.failures.length > 0) {
    for (var resultIndex = 0; resultIndex < result.failures.length; resultIndex += 1) {
        process.stderr.write("FAIL: " + result.failures[resultIndex] + "\n");
    }
    process.exitCode = 1;
} else {
    process.stdout.write(
        "Code-generator core checks passed: " + result.checks + " assertions.\n"
    );
}

if (process.argv[2]) {
    var outputDirectory = path.resolve(process.argv[2]);
    fs.mkdirSync(outputDirectory, { recursive: true });
    fs.writeFileSync(
        path.join(outputDirectory, "generated_smoke.jl"),
        result.bundle,
        "utf8"
    );
    process.stdout.write(
        "Wrote generated Julia smoke bundle to " +
            path.join(outputDirectory, "generated_smoke.jl") +
            "\n"
    );
}

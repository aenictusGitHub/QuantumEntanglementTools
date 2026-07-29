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
for (var idIndex = 0; idIndex < result.requiredDomIds.length; idIndex += 1) {
    var id = result.requiredDomIds[idIndex];
    if (pageSource.indexOf('id="' + id + '"') === -1) {
        result.failures.push("documentation page is missing required DOM id " + id);
    }
}
var referencedIdPattern = /byId\("([^"]+)"\)/g;
var referencedIdMatch;
while ((referencedIdMatch = referencedIdPattern.exec(uiSource)) !== null) {
    if (pageSource.indexOf('id="' + referencedIdMatch[1] + '"') === -1) {
        result.failures.push(
            "UI adapter references missing DOM id " + referencedIdMatch[1]
        );
    }
}

if (result.failures.length > 0) {
    for (var failureIndex = 0; failureIndex < result.failures.length; failureIndex += 1) {
        process.stderr.write("FAIL: " + result.failures[failureIndex] + "\n");
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

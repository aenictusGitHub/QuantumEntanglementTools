# Optional extension template

This directory is a design placeholder, not a loadable extension.

Before adding an adapter, follow
`docs/src/integrating_external_packages.md` and ADR 0004. A concrete extension
contribution should add:

```text
ext/QuantumEntanglementToolsBackendNameExt.jl
test/extensions/backend_name/
docs/src/backends/backend_name.md
```

It must also update the root project's `[weakdeps]`, `[extensions]`, and
`[compat]` entries; report capabilities and backend version; translate only
documented public backend APIs into package-owned result types; and test both
load orders plus absence of the dependency.

Do not copy this placeholder into runtime code without replacing every abstract
name, auditing global side effects, and adding a dedicated test environment.

For a reviewed implementation, see the EntanglementDetection adapter in
`ext/`, its dedicated environment under
`test/extensions/entanglement_detection/`, and
`docs/src/entanglement_detection_extension.md`.

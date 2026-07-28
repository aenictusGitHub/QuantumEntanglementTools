# Integrating external packages

External packages integrate through Julia weak dependencies and package
extensions. The core API and result semantics remain owned by
`QuantumEntanglementTools`.

## Adapter checklist

1. Identify a concrete user operation and public backend API.
2. Audit version, license, maintenance, return values, exceptions, RNG, logging,
   stdout, threading, and global state.
3. Document basis order, dimensions, normalization, allocation, and certificate
   differences.
4. Add `[weakdeps]` and `[extensions]` metadata without adding a core runtime
   dependency.
5. Implement only public methods and contain raw backend types.
6. Test absence of the backend, both load orders, repeated loading,
   precompilation, failures, state preservation, and method ambiguities.
7. Publish backend capabilities and version metadata.

An illustrative extension layout is:

```text
ext/
  QuantumEntanglementToolsExampleBackendExt.jl
test/
  extensions/
    example_backend/
```

The extension module imports both packages and adds narrowly scoped methods to
existing core generic functions. It must not re-export the backend namespace or
install methods dynamically. The contributor-facing placeholder is in
`dev/extension_template/README.md`.

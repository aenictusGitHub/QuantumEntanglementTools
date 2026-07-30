# MATLAB compatibility API

These entry points preserve only the reviewed QETLAB call shapes recorded in
the migration ledger. They do not imply complete QETLAB or general MATLAB
equivalence. In particular, `AbsPPTConstraints` preserves positional
`DIM`, `ESC_IF_NPOS`, and `LIM`, while `IsAbsPPT` preserves the upstream
`1/0/-1` tri-state surface without collapsing native numerical or resource
boundaries.

```@autodocs
Modules = [QuantumEntanglementTools.MATLABCompat]
Private = false
```

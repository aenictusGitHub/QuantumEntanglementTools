# Unextendible product-basis catalog

`upb` constructs every executable family in the pinned QETLAB `UPB` catalog.
It returns one stable `UPBConstruction` object instead of changing the output
shape according to the number of requested MATLAB outputs.

```julia
construction = upb(:tiles)

construction.family       # :tiles
construction.dimensions   # (3, 3)
construction.cardinality  # 5
construction.local_factors[1][:, 2]
construction.global_vectors[:, 2]
```

For state $j$, the global column is

```math
u_j =
u_j^{(1)} \otimes u_j^{(2)} \otimes \cdots \otimes u_j^{(p)}.
```

Subsystem 1 is the slowest-varying tensor factor. The local factors and global
columns are always both present, and `dimensions` remains in the caller's
original subsystem order.

## Result and diagnostics

The result records:

- local factors, global vectors, dimensions, cardinality, and canonical family;
- the primary reference key, citation, and URL;
- whether the construction is closed-form, randomized full-spark, or a
  complete product basis;
- normalization, pairwise-orthogonality, and tensor-reconstruction residuals;
- explicit RNG use, attempts, random draws, work, and full-spark minors;
- the entry, work, minor, and retry limits used for the call;
- the arithmetic precision and the disposition of symbolic simplification.

`verification_kind == :tolerance_robust` is a numerical construction check. It
does not relabel floating entries as exact arithmetic. Computational full bases
use `:exact_structure` because their zeros and ones establish the relevant
structure exactly in each supported floating type.

The constructor does not silently normalize caller data: it has no caller
vectors to repair. Family formulas that mathematically specify normalized
vectors are evaluated as written, then independently checked. A result is not
returned if normalization or orthogonality falls in or beyond the configured
boundary band.

## Named-family coverage

The source names are case-insensitive. Julia-native underscored symbols and the
pinned spellings are both accepted.

| Pinned name | Native family | Parameters | Dimensions | Cardinality | Construction |
|---|---|---:|---:|---:|---|
| `Pyramid` | `:pyramid` | none | `(3,3)` | 5 | closed form |
| `Tiles` | `:tiles` | none | `(3,3)` | 5 | closed form |
| `GenTiles1` | `:generalized_tiles_1` | even `d ≥ 4` | `(d,d)` | `(d-1)^2` | closed form |
| `GenTiles2` | `:generalized_tiles_2` | `m,n`, `n ≥ m ≥ 3`, `n ≥ 4` | `(m,n)` | `mn-2m+1` | closed form |
| `Min4x4` | `:minimum_4x4` | none | `(4,4)` | 8 | closed form |
| `QuadRes` | `:quad_residue` | odd `d`, `2d-1` prime | `(d,d)` | `2d-1` | closed form |
| `SixParam` | `:six_parameter` | six angles | `(3,3)` | 5 | closed form |
| `Shifts` | `:shifts` | none | `(2,2,2)` | 4 | closed form |
| `GenShifts` | `:generalized_shifts` | odd `p ≥ 3` | `p` qubits | `p+1` | closed form |
| `Feng2x2x3` | `:feng_2x2x3` | none | `(2,2,3)` | 6 | closed form |
| `Feng2x2x5` | `:feng_2x2x5` | none | `(2,2,5)` | 8 | closed form |
| `Feng2x2x2x2` | `:feng_2x2x2x2` | none | four qubits | 6 | closed form |
| `Feng4x4` | `:feng_4x4` | none | `(4,4)` | 8 | closed form |
| `Feng2x2x2x4` | `:feng_2x2x2x4` | none | `(2,2,2,4)` | 8 | closed form |
| `Feng2x2x2x2x5` | `:feng_2x2x2x2x5` | none | `(2,2,2,2,5)` | 10 | closed form |
| `John2^8` | `:johnston_2_power_8` | none | eight qubits | 11 | closed form |
| `John2^4k` | `:johnston_2_power_4k` | `p ≥ 8`, `p = 0 mod 4` | `p` qubits | `p+4` | corrected closed form |
| `CJBip46` | `:chen_johnston_4x6` | none | `(4,6)` | 10 | closed form |
| `AlonLovasz` | `:alon_lovasz` | dimensions | routed dimensions | counting lower bound | randomized full spark |
| `CJBip` | `:chen_johnston_bipartite` | dimensions | routed dimensions | `2d_max` | randomized full spark |
| `CJ4k1` | `:chen_johnston_4k1` | `d = 1 mod 4`, `d ≥ 5` | `(2,2,d)` | `d+3` | randomized full spark |

The last three names expose safe expert access to construction branches that
the pinned file normally reaches through numeric dimension dispatch.

## Dimension dispatch

Passing dimensions follows the pinned ordering and then restores the input
subsystem order:

```julia
construction = upb((4, 2, 2, 2))

construction.family
# :feng_2x2x2x4

construction.dimensions
# (4, 2, 2, 2)
```

The routed branches are:

1. a complete basis for one party or a bipartite space with a local dimension
   at most two;
2. Tiles, Feng `4×4`, and quadratic-residue equal-dimensional families;
3. the fixed three-, four-, five-, and eight-party families in the table;
4. Johnston's `4k` qubit family and generalized Shifts for odd qubit counts;
5. Alon–Lovász, Chen–Johnston largest-dimension equality, fixed `4×6`, and
   Chen–Johnston `(2,2,4k+1)` constructions;
6. a structured `UPBConstructionUnavailable` exception otherwise.

The exception distinguishes `:minimum_unknown` from
`:known_but_not_in_catalog`. The latter retains the exact theorem-backed
minimum. This preserves the public catalog's intentional hard-construction
boundary without inventing vectors:

```julia
try
    upb((2, 2, 2, 2, 2, 2))
catch err
    @assert err isa UPBConstructionUnavailable
    @assert err.reason == :known_but_not_in_catalog
    @assert err.known_minimum == 8
end
```

## Explicit RNG and bounded searches

Randomized routes require an explicit `AbstractRNG`:

```julia
using Random

rng = MersenneTwister(2026)
construction = upb(rng, (3, 4))

construction.family          # :alon_lovasz
construction.rng_used        # true
construction.attempts
construction.minors_checked
```

There is no implicit default-RNG fallback. Supplying the same RNG type and seed
reproduces the same local factors. Library code never seeds or reads the global
random stream.

The upstream retry loops and minor enumeration have no general resource
bound. The Julia API instead exposes:

```julia
upb(
    rng,
    (6, 6);
    max_attempts=64,
    max_minors=250_000,
    max_work=50_000_000,
    max_local_entries=1_000_000,
    max_global_entries=2_000_000,
)
```

Entry guards run before the corresponding local or global allocation.
Full-spark minor counts are checked before enumeration. A deterministic limit
raises `UPBResourceLimitError` with the resource, required amount when known,
and configured limit. A bounded random search that exhausts all attempts raises
`UPBConstructionUnavailable` with reason `:randomized_search_exhausted`.

Randomized rank decisions currently use LAPACK SVD and therefore accept
`Float32` and `Float64`. Deterministic families additionally accept
`BigFloat`:

```julia
setprecision(BigFloat, 192) do
    high_precision = upb(:quad_residue, 3; real_type=BigFloat)
    @assert eltype(high_precision.global_vectors) == Complex{BigFloat}
end
```

## Independent validation

The focused property suite validates every family for:

- local and global dimensions and expected cardinality;
- finite normalized local vectors;
- columnwise tensor-product reconstruction;
- pairwise global orthogonality;
- unextendibility via the independent bounded `is_upb` partition/rank search
  whenever exhaustive enumeration is practical;
- the independent qubit ray-cover graph criterion for the two large Johnston
  catalogs;
- RNG reproducibility, global-stream isolation, generic precision, invalid
  inputs, and every resource guard.

The round-robin private `_upbc_one_factorization` replaces the pinned helper
without capability loss. Tests prove that for 2, 4, 6, and 8 labelled vertices
every unordered edge appears exactly once and every row is a perfect matching.
The helper remains private because it exists only to support the public
Johnston construction.

The committed source-free Octave 11.3 fixture compares 17 deterministic
families with the pinned QETLAB revision up to independent local-vector phases.
Randomized families use theorem, full-spark, orthogonality, and independent
unextendibility properties instead of pretending that different random draws
must agree entrywise.

## Reviewed upstream defects and deviations

### `GenTiles1(2)`

The pinned condition rejects odd dimensions but accepts `2`. Its loops then
return only the all-ones product vector, which is extendible. The native domain
requires an even dimension of at least four, matching the mathematical family.

### `John2^4k` reshape order

The pinned final one-factorization branch reshapes a `2×2×(2k+2)` basis array
in column-major order and then indexes it as though all first basis vectors
preceded all orthogonal partners. For eight parties, the pinned product Gram
matrix has off-diagonal residual approximately
`0.13778053319420595`, so those returned vectors are not pairwise orthogonal.
The committed oracle retains this minimized behavior.

The native branch constructs the orthogonality graphs in Lemma 4 of Johnston's
paper directly: the $B_{j,k}$ complete-bipartite components are followed by a
round-robin one-factorization of the two complete-graph components. It returns
an orthogonal, independently unextendible `p+4` family for every supported
`p`, including the nonminimal 12-state eight-qubit member.

### Symbolic `NICE`

The Symbolic Math Toolbox code belongs to nested Chen–Johnston helper
functions. The top-level pinned `UPB` routes call those helpers without
forwarding `NICE`, so the symbolic branch is unreachable through the executable
public function. The native numeric construction therefore loses no public
capability and records
`symbolic_simplification == :not_required_for_executable_public_route`.

## MATLAB compatibility

The compatibility wrapper has an explicit output choice because Julia has no
MATLAB `nargout` dispatch:

```julia
MATLABCompat.UPB(
    input,
    args...;
    output=:global,      # :global, :local, or :structured
    kwargs...,
)

MATLABCompat.UPB(
    rng,
    input,
    args...;
    output=:structured,
    kwargs...,
)
```

- `output=:global` returns a copy of `global_vectors`;
- `output=:local` returns the tuple of local-factor matrices;
- `output=:structured` returns the full `UPBConstruction`;
- randomized routes require the leading explicit RNG;
- numeric `DIM` calls preserve the optional positional `VERBOSE` reference
  flag;
- references are always retained in structured output and `VERBOSE=true`
  additionally prints the selected citation;
- symbolic `NICE` is documented as unreachable rather than accepted and
  ignored.

The compatibility wrapper should preserve pinned family spelling and numeric
dimension routing, but must retain the two mathematical corrections above.

## Provenance and benchmark plan

The row-level provenance entry should identify:

- pinned `UPB.m` SHA-256
  `d60d9d2375d773c7e59930991ddf73427c9ffb0c77ded640b7f635f99ec13544`;
- pinned `helpers/one_factorization.m` SHA-256
  `7b042b3729b96dc4646c57adc24dcd5e4f59cf710520b0a6f98e7f3fd0c531ca`;
- source-informed independent implementation, the primary sources below,
  focused property tests, and the phase-aware oracle fixture;
- the two corrected source defects and the unreachable symbolic-helper
  disposition.

Quick benchmark coverage should include:

- closed-form Tiles;
- generalized Tiles 1 in dimension 8;
- Johnston's 12-qubit construction;
- seeded Alon–Lovász `(3,4)`;
- seeded Chen–Johnston `(6,6)`;
- local-factor construction separately from required global-vector
  materialization where the harness supports setup separation.

Benchmarks are local smoke evidence only until a versioned regression artifact
records machine, Julia, BLAS, input, runtime, allocations, and peak memory.

## Primary references

- [Bennett et al., *Unextendible Product Bases and Bound
  Entanglement*](https://arxiv.org/abs/quant-ph/9808030).
- [DiVincenzo et al., *Unextendible Product Bases, Uncompletable Product
  Bases and Bound Entanglement*](https://arxiv.org/abs/quant-ph/9908070).
- [Alon and Lovász, *Unextendible Product
  Bases*](https://doi.org/10.1006/jcta.2000.3122).
- [Feng, *Unextendible Product Bases and 1-Factorization of Complete
  Graphs*](https://doi.org/10.1016/j.dam.2005.10.011).
- [Johnston, *The Minimum Size of Qubit Unextendible Product
  Bases*](https://arxiv.org/abs/1302.1604).
- [Chen and Johnston, *The Minimum Size of Unextendible Product Bases in
  the Bipartite Case (and Some Multipartite
  Cases)*](https://arxiv.org/abs/1301.1406).

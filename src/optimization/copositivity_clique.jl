# Source-informed independent Julia implementations based on the executable
# contracts of QETLAB IsCopositive.m and CliqueNumber.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
#
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# Numerical hierarchy values are deliberately kept separate from exact
# certificates. In particular, a floating SOS lower bound close to zero is
# never promoted to a proof of copositivity, and a transformed floating
# Motzkin--Straus value is never promoted to a certified clique upper bound.

using LinearAlgebra
using Random: AbstractRNG
using SparseArrays

"""
    CopositivityStatus

Certificate-aware outcome of [`copositivity_criterion`](@ref). Only
`CopositivityCertifiedTrue` and `CopositivityCertifiedFalse` carry Boolean
verdicts. Solver bounds, unavailable backends, limits, and tolerance
boundaries remain inconclusive.
"""
@enum CopositivityStatus::UInt8 begin
    CopositivityCertifiedTrue
    CopositivityCertifiedFalse
    CopositivityHierarchyUnknown
    CopositivityNumericalBoundary
    CopositivityBackendUnavailable
    CopositivityResourceLimit
    CopositivityBackendFailure
end

"""
    CopositivityWitness

Exact nonnegative-simplex witness for failure of copositivity.
`simplex_vector` and `value` use `Rational{BigInt}` arithmetic, and
`value == simplex_vector' * matrix * simplex_vector < 0` for the matrix
entries exactly as supplied. `source` records whether the witness was a
coordinate ray, a two-coordinate ray, the uniform ray, or a rationalized
sample. The coordinate storage is owned by the result and read-only through
the public array interface. Instances are produced by
[`copositivity_criterion`](@ref), not by a public unchecked constructor.
"""
struct CopositivityWitness{V,Q<:Rational}
    simplex_vector::V
    value::Q
    source::Symbol

    function CopositivityWitness(
        token::_ValidatedConstructorToken, simplex_vector::V, value::Q, source::Symbol
    ) where {V,Q<:Rational}
        _require_validated_constructor_token(token)
        simplex_vector isa _ReadOnlyPlanArray || throw(
            ArgumentError(
                "copositivity witness coordinates must use owned read-only storage"
            ),
        )
        eltype(simplex_vector) === Rational{BigInt} || throw(
            ArgumentError("copositivity witness coordinates must use Rational{BigInt}")
        )
        !isempty(simplex_vector) ||
            throw(ArgumentError("a copositivity witness must be nonempty"))
        all(coordinate -> coordinate >= 0, simplex_vector) ||
            throw(ArgumentError("a copositivity witness must be entrywise nonnegative"))
        sum(simplex_vector) == 1 ||
            throw(ArgumentError("a copositivity witness must lie on the simplex"))
        value < 0 || throw(ArgumentError("a copositivity witness value must be negative"))
        source in
        (:coordinate_ray, :two_coordinate_ray, :uniform_ray, :rationalized_sample) ||
            throw(ArgumentError("unsupported copositivity witness source $source"))
        return new{V,Q}(simplex_vector, value, source)
    end
end

"""
    CopositivityResult

Status-rich result for the minimum of `y' * C * y` over the nonnegative
simplex. `lower_bound` is supplied by the selected hierarchy and
`upper_bound` is attained by an explicit point. Floating hierarchy bounds
are numerical evidence and do not by themselves set `verdict`.

`verdict` is present only for an exact entrywise-nonnegative or exact-PSD
sufficient certificate, or for an exact rational witness with negative
value outside the requested tolerance band. `polynomial === nothing` when an
exact branch decides the result before hierarchy construction. Instances are
produced by [`copositivity_criterion`](@ref), not by a public unchecked
constructor.
"""
struct CopositivityResult{L,U,W,D,P,H,T}
    status::CopositivityStatus
    verdict::Union{Nothing,Bool}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    lower_bound::L
    upper_bound::U
    lower_kind::Symbol
    upper_kind::Symbol
    witness::W
    psd_diagnostic::D
    polynomial::P
    hierarchy::Symbol
    hierarchy_level::Int
    hierarchy_result::H
    samples_requested::Int
    samples_evaluated::Int
    tolerance::T
    message::String

    function CopositivityResult(
        token::_ValidatedConstructorToken,
        status::CopositivityStatus,
        verdict::Union{Nothing,Bool},
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        lower_bound,
        upper_bound,
        lower_kind::Symbol,
        upper_kind::Symbol,
        witness,
        psd_diagnostic,
        polynomial,
        hierarchy::Symbol,
        hierarchy_level::Int,
        hierarchy_result,
        samples_requested::Int,
        samples_evaluated::Int,
        tolerance,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        hierarchy in (:sos, :nosdp) ||
            throw(ArgumentError("unsupported copositivity hierarchy $hierarchy"))
        hierarchy_level >= 0 || throw(ArgumentError("hierarchy_level must be nonnegative"))
        0 <= samples_evaluated <= samples_requested ||
            throw(ArgumentError("sample counters are inconsistent"))
        tolerance isa Real &&
        !(tolerance isa Bool) &&
        isfinite(tolerance) &&
        tolerance >= 0 ||
            throw(ArgumentError("tolerance must be a finite nonnegative real number"))
        for (value, name) in ((lower_bound, "lower_bound"), (upper_bound, "upper_bound"))
            value === nothing && continue
            value isa Real && !(value isa Bool) && isfinite(value) ||
                throw(ArgumentError("$name must be finite when present"))
        end

        if status === CopositivityCertifiedTrue
            verdict === true && certified && certificate_kind !== nothing || throw(
                ArgumentError(
                    "CopositivityCertifiedTrue requires a certified true verdict and certificate",
                ),
            )
            witness === nothing ||
                throw(ArgumentError("a positive copositivity certificate has no witness"))
        elseif status === CopositivityCertifiedFalse
            verdict === false &&
            certified &&
            certificate_kind !== nothing &&
            witness isa CopositivityWitness || throw(
                ArgumentError(
                    "CopositivityCertifiedFalse requires a certified false verdict and witness",
                ),
            )
            upper_bound == witness.value || throw(
                ArgumentError(
                    "a negative certificate upper bound must equal its witness value"
                ),
            )
        else
            verdict === nothing && !certified && certificate_kind === nothing || throw(
                ArgumentError(
                    "inconclusive copositivity statuses cannot carry a verdict or certificate",
                ),
            )
            witness === nothing || throw(
                ArgumentError("an inconclusive copositivity result cannot carry a witness"),
            )
        end
        return new{
            typeof(lower_bound),
            typeof(upper_bound),
            typeof(witness),
            typeof(psd_diagnostic),
            typeof(polynomial),
            typeof(hierarchy_result),
            typeof(tolerance),
        }(
            status,
            verdict,
            certified,
            certificate_kind,
            lower_bound,
            upper_bound,
            lower_kind,
            upper_kind,
            witness,
            psd_diagnostic,
            polynomial,
            hierarchy,
            hierarchy_level,
            hierarchy_result,
            samples_requested,
            samples_evaluated,
            tolerance,
            String(message),
        )
    end
end

function Base.show(io::IO, result::CopositivityResult)
    return print(
        io,
        "CopositivityResult(status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", lower_bound=",
        result.lower_bound,
        ", upper_bound=",
        result.upper_bound,
        ")",
    )
end

"""
    CliqueNumberStatus

Status of the polynomial work retained by [`clique_number_bounds`](@ref).
The returned discrete interval is always backed by exact graph-theoretic
certificates, even when the optional hierarchy is unavailable or fails.
"""
@enum CliqueNumberStatus::UInt8 begin
    CliqueNumberExact
    CliqueNumberNumericalHierarchy
    CliqueNumberBackendUnavailable
    CliqueNumberResourceLimit
    CliqueNumberBackendFailure
end

"""
    MotzkinStrausWitness

An exact rational point on the nonnegative simplex. The Motzkin--Straus
theorem turns `value` into the certified integer `implied_lower_bound`.
The simplex storage is owned by the result and read-only through the public
array interface. Instances are produced by [`clique_number_bounds`](@ref),
not by a public unchecked constructor.
"""
struct MotzkinStrausWitness{V,Q<:Rational}
    simplex_vector::V
    value::Q
    implied_lower_bound::Int

    function MotzkinStrausWitness(
        token::_ValidatedConstructorToken,
        simplex_vector::V,
        value::Q,
        implied_lower_bound::Int,
    ) where {V,Q<:Rational}
        _require_validated_constructor_token(token)
        simplex_vector isa _ReadOnlyPlanArray || throw(
            ArgumentError(
                "Motzkin--Straus witness coordinates must use owned read-only storage"
            ),
        )
        eltype(simplex_vector) === Rational{BigInt} || throw(
            ArgumentError("Motzkin--Straus witness coordinates must use Rational{BigInt}"),
        )
        !isempty(simplex_vector) ||
            throw(ArgumentError("a Motzkin--Straus witness must be nonempty"))
        all(coordinate -> coordinate >= 0, simplex_vector) ||
            throw(ArgumentError("a Motzkin--Straus witness must be entrywise nonnegative"))
        sum(simplex_vector) == 1 ||
            throw(ArgumentError("a Motzkin--Straus witness must lie on the simplex"))
        0 <= value < 1 ||
            throw(ArgumentError("a Motzkin--Straus witness value must lie in [0,1)"))
        reciprocal = inv(1 - value)
        expected = cld(numerator(reciprocal), denominator(reciprocal))
        implied_lower_bound == expected || throw(
            ArgumentError(
                "the Motzkin--Straus implied lower bound is inconsistent with its value"
            ),
        )
        implied_lower_bound >= 1 ||
            throw(ArgumentError("a clique lower bound must be positive"))
        return new{V,Q}(simplex_vector, value, implied_lower_bound)
    end
end

"""
    CliqueNumberResult

Certified integer bounds on a simple undirected graph's clique number.
`best_clique` proves a combinatorial lower bound. The edge count, maximum
degree, and deterministic greedy coloring prove the structural upper bound.
An exact rationalized Motzkin--Straus sample may improve the lower bound.

`upper_certificate` retains every component upper bound, the maximum degree,
and the full proper coloring. `upper_certificate_kinds` lists exactly the
components attaining `upper_bound`. The clique, coloring, and any
Motzkin--Straus simplex use owned storage that is read-only through the public
array interface.

`uncertified_upper_candidate` is the integer suggested by a floating
hierarchy outer value after outward tolerance padding. It is intentionally
not used to tighten `upper_bound`. `polynomial === nothing` when the
graph-theoretic certificates are already exact. Instances are produced by
[`clique_number_bounds`](@ref), not by a public unchecked constructor.
"""
struct CliqueNumberResult{P,H,CL,CU,N,W,C,T}
    status::CliqueNumberStatus
    lower_bound::Int
    upper_bound::Int
    exact::Bool
    bounds_certified::Bool
    lower_certificate_kind::Symbol
    upper_certificate_kinds::Tuple{Vararg{Symbol}}
    upper_certificate::C
    best_clique::_ReadOnlyPlanVector
    vertex_count::Int
    edge_count::BigInt
    polynomial::P
    hierarchy::Symbol
    hierarchy_level::Int
    hierarchy_result::H
    continuous_lower_bound::CL
    continuous_upper_bound::CU
    uncertified_upper_candidate::N
    motzkin_straus_witness::W
    samples_requested::Int
    samples_evaluated::Int
    tolerance::T
    message::String

    function CliqueNumberResult(
        token::_ValidatedConstructorToken,
        status::CliqueNumberStatus,
        lower_bound::Int,
        upper_bound::Int,
        exact::Bool,
        bounds_certified::Bool,
        lower_certificate_kind::Symbol,
        upper_certificate_kinds::Tuple{Vararg{Symbol}},
        upper_certificate,
        best_clique::_ReadOnlyPlanVector,
        vertex_count::Int,
        edge_count::BigInt,
        polynomial,
        hierarchy::Symbol,
        hierarchy_level::Int,
        hierarchy_result,
        continuous_lower_bound,
        continuous_upper_bound,
        uncertified_upper_candidate,
        motzkin_straus_witness,
        samples_requested::Int,
        samples_evaluated::Int,
        tolerance,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        vertex_count >= 1 || throw(ArgumentError("vertex_count must be positive"))
        1 <= lower_bound <= upper_bound <= vertex_count ||
            throw(ArgumentError("clique bounds must lie in 1:vertex_count"))
        exact == (lower_bound == upper_bound) ||
            throw(ArgumentError("exact must agree with equality of the certified bounds"))
        (status === CliqueNumberExact) == exact ||
            throw(ArgumentError("CliqueNumberExact must agree with the exact flag"))
        bounds_certified ||
            throw(ArgumentError("CliqueNumberResult bounds must be graph-certified"))
        hierarchy in (:sos, :nosdp) ||
            throw(ArgumentError("unsupported clique hierarchy $hierarchy"))
        hierarchy_level >= 0 || throw(ArgumentError("hierarchy_level must be nonnegative"))
        0 <= samples_evaluated <= samples_requested ||
            throw(ArgumentError("sample counters are inconsistent"))
        tolerance isa Real &&
        !(tolerance isa Bool) &&
        isfinite(tolerance) &&
        tolerance >= 0 ||
            throw(ArgumentError("tolerance must be a finite nonnegative real number"))
        edge_count >= 0 || throw(ArgumentError("edge_count must be nonnegative"))
        edge_count <= BigInt(vertex_count) * (vertex_count - 1) ÷ 2 ||
            throw(ArgumentError("edge_count exceeds the simple-graph maximum"))
        !isempty(best_clique) ||
            throw(ArgumentError("best_clique must contain at least one vertex"))
        issorted(best_clique) && allunique(best_clique) ||
            throw(ArgumentError("best_clique must contain sorted distinct vertex indices"))
        all(vertex -> 1 <= vertex <= vertex_count, best_clique) ||
            throw(ArgumentError("best_clique contains an out-of-range vertex"))
        length(best_clique) <= lower_bound ||
            throw(ArgumentError("best_clique is longer than the claimed lower bound"))
        lower_certificate_kind in (:greedy_clique, :motzkin_straus_witness) ||
            throw(ArgumentError("unsupported clique lower-certificate kind"))
        if lower_certificate_kind === :greedy_clique
            length(best_clique) == lower_bound || throw(
                ArgumentError("a greedy-clique lower bound must equal best_clique length"),
            )
        else
            motzkin_straus_witness isa MotzkinStrausWitness || throw(
                ArgumentError("a Motzkin--Straus lower bound requires its exact witness"),
            )
            motzkin_straus_witness.implied_lower_bound == lower_bound || throw(
                ArgumentError("the Motzkin--Straus witness does not imply lower_bound")
            )
        end

        required_keys = (
            :edge_count_bound,
            :maximum_degree,
            :maximum_degree_bound,
            :greedy_coloring_bound,
            :coloring,
            :attaining_kinds,
        )
        all(key -> hasproperty(upper_certificate, key), required_keys) ||
            throw(ArgumentError("the clique upper certificate is incomplete"))
        coloring = upper_certificate.coloring
        coloring isa _ReadOnlyPlanVector ||
            throw(ArgumentError("the clique coloring must use owned read-only storage"))
        length(coloring) == vertex_count ||
            throw(ArgumentError("the clique coloring has the wrong vertex count"))
        all(color -> color >= 1, coloring) ||
            throw(ArgumentError("clique coloring labels must be positive"))
        maximum(coloring) == upper_certificate.greedy_coloring_bound || throw(
            ArgumentError("the clique coloring does not attain its recorded color bound"),
        )
        upper_certificate.edge_count_bound ==
        _copclique_edge_upper(edge_count, vertex_count) ||
            throw(ArgumentError("the edge-count upper bound is inconsistent"))
        0 <= upper_certificate.maximum_degree < vertex_count ||
            throw(ArgumentError("maximum_degree is outside the simple-graph range"))
        upper_certificate.maximum_degree_bound ==
        min(vertex_count, upper_certificate.maximum_degree + 1) ||
            throw(ArgumentError("the maximum-degree upper bound is inconsistent"))
        upper_bound == min(
            upper_certificate.edge_count_bound,
            upper_certificate.maximum_degree_bound,
            upper_certificate.greedy_coloring_bound,
        ) || throw(ArgumentError("upper_bound does not match its component certificates"))
        upper_certificate_kinds == upper_certificate.attaining_kinds ||
            throw(ArgumentError("upper_certificate_kinds is inconsistent"))
        isempty(upper_certificate_kinds) &&
            throw(ArgumentError("at least one upper certificate must attain upper_bound"))
        for (value, name) in (
            (continuous_lower_bound, "continuous_lower_bound"),
            (continuous_upper_bound, "continuous_upper_bound"),
        )
            value === nothing && continue
            value isa Real && !(value isa Bool) && isfinite(value) ||
                throw(ArgumentError("$name must be finite when present"))
        end
        if uncertified_upper_candidate !== nothing
            uncertified_upper_candidate isa Int &&
            1 <= uncertified_upper_candidate <= vertex_count || throw(
                ArgumentError("uncertified_upper_candidate is outside the graph range")
            )
        end
        return new{
            typeof(polynomial),
            typeof(hierarchy_result),
            typeof(continuous_lower_bound),
            typeof(continuous_upper_bound),
            typeof(uncertified_upper_candidate),
            typeof(motzkin_straus_witness),
            typeof(upper_certificate),
            typeof(tolerance),
        }(
            status,
            lower_bound,
            upper_bound,
            exact,
            bounds_certified,
            lower_certificate_kind,
            upper_certificate_kinds,
            upper_certificate,
            best_clique,
            vertex_count,
            edge_count,
            polynomial,
            hierarchy,
            hierarchy_level,
            hierarchy_result,
            continuous_lower_bound,
            continuous_upper_bound,
            uncertified_upper_candidate,
            motzkin_straus_witness,
            samples_requested,
            samples_evaluated,
            tolerance,
            String(message),
        )
    end
end

function Base.show(io::IO, result::CliqueNumberResult)
    return print(
        io,
        "CliqueNumberResult(status=",
        result.status,
        ", bounds=",
        result.lower_bound,
        ":",
        result.upper_bound,
        ", exact=",
        result.exact,
        ")",
    )
end

const _COPCLIQUE_DEFAULT_MAX_MATRIX_DIMENSION = 128
const _COPCLIQUE_DEFAULT_MAX_MATRIX_ENTRIES = 16_384
const _COPCLIQUE_DEFAULT_MAX_GREEDY_WORK = 5_000_000
const _COPCLIQUE_DEFAULT_MAX_CERTIFICATE_BITS = 100_000

function _copclique_nonnegative_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _copclique_positive_integer(value, name::AbstractString)
    checked = _copclique_nonnegative_integer(value, name)
    checked > 0 || throw(ArgumentError("$name must be a positive integer"))
    return checked
end

function _copclique_optional_limit(value, name::AbstractString; allow_zero::Bool=false)
    value === nothing && return nothing
    return BigInt(
        if allow_zero
            _copclique_nonnegative_integer(value, name)
        else
            _copclique_positive_integer(value, name)
        end,
    )
end

function _copclique_hierarchy(value)
    value isa Symbol ||
        throw(ArgumentError("hierarchy must be :sos or :nosdp; got $(repr(value))"))
    value in (:sos, :nosdp) ||
        throw(ArgumentError("hierarchy must be :sos or :nosdp; got $(repr(value))"))
    return value
end

function _copclique_validate_backend(
    hierarchy::Symbol, backend::AbstractOptimizationBackend
)
    hierarchy === :nosdp &&
        !(backend isa NoOptimizationBackend) &&
        throw(
            ArgumentError(
                "the :nosdp hierarchy is solver-free; do not pass an optimization backend"
            ),
        )
    return backend
end

function _copclique_matrix_contract(
    matrix::AbstractMatrix; graph::Bool, max_matrix_dimension, max_matrix_entries
)
    Base.require_one_based_indexing(matrix)
    rows, columns = size(matrix)
    rows == columns || throw(
        DimensionMismatch(
            if graph
                "an adjacency matrix must be square; got $(size(matrix))"
            else
                "a copositivity matrix must be square; got $(size(matrix))"
            end,
        ),
    )
    rows > 0 || throw(ArgumentError(
        if graph
            "an adjacency matrix must contain at least one vertex"
        else
            "a copositivity matrix must be nonempty"
        end,
    ))
    dimension_limit = _copclique_optional_limit(
        max_matrix_dimension, "max_matrix_dimension"
    )
    dimension_limit !== nothing &&
        rows > dimension_limit &&
        throw(
            ArgumentError(
                "matrix dimension $rows exceeds max_matrix_dimension=$dimension_limit"
            ),
        )
    entries = BigInt(rows)^2
    entry_limit = _copclique_optional_limit(
        max_matrix_entries, "max_matrix_entries"; allow_zero=true
    )
    entry_limit !== nothing &&
        entries > entry_limit &&
        throw(
            ArgumentError(
                "matrix validation needs $entries logical entries, exceeding " *
                "max_matrix_entries=$entry_limit",
            ),
        )

    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Real || throw(
        ArgumentError(
            if graph
                "an adjacency matrix needs a concrete real element type; got $value_type"
            else
                "a copositivity matrix needs a concrete real element type; got $value_type"
            end,
        ),
    )
    supported_type =
        value_type <: Integer || value_type <: Rational || value_type <: AbstractFloat
    supported_type || throw(
        ArgumentError(
            "matrix entries must use integer, rational, or floating-point storage; " *
            "got $value_type",
        ),
    )
    !graph &&
        value_type === Bool &&
        throw(ArgumentError("Bool is not a supported copositivity coefficient type"))

    maximum_magnitude = zero(value_type)
    for column in 1:columns, row in 1:rows
        value = matrix[row, column]
        isfinite(value) || throw(ArgumentError(
            if graph
                "adjacency entry ($row,$column) must be finite"
            else
                "matrix entry ($row,$column) must be finite"
            end,
        ))
        if graph
            (iszero(value) || value == one(value)) || throw(
                ArgumentError(
                    "adjacency entry ($row,$column) is $(repr(value)); " *
                    "only exact zero/one entries are accepted",
                ),
            )
            row == column &&
                !iszero(value) &&
                throw(
                    ArgumentError("self-loops are not accepted; diagonal entry $row is one")
                )
        else
            maximum_magnitude = max(maximum_magnitude, abs(value))
        end
    end
    issymmetric(matrix) || throw(
        ArgumentError(
            if graph
                "the adjacency matrix must be exactly symmetric; no repair is applied"
            else
                "the matrix must be exactly symmetric; no Hermitian part is substituted"
            end,
        ),
    )
    exact = value_type <: Integer || value_type <: Rational
    return rows, value_type, exact, maximum_magnitude
end

function _copclique_tolerance(value_type, exact::Bool, maximum_magnitude; atol, rtol)
    float_type = typeof(float(zero(value_type)))
    float_type <: AbstractFloat ||
        throw(ArgumentError("cannot choose a real floating tolerance type for $value_type"))
    absolute = atol === nothing ? zero(float_type) : atol
    relative = rtol === nothing ? (exact ? zero(float_type) : sqrt(eps(float_type))) : rtol
    for (value, name) in ((absolute, "atol"), (relative, "rtol"))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    converted_absolute = convert(float_type, absolute)
    converted_relative = convert(float_type, relative)
    exact &&
        (!iszero(converted_absolute) || !iszero(converted_relative)) &&
        throw(
            ArgumentError(
                "exact integer and rational inputs use exact comparisons; " *
                "atol and rtol must both be zero",
            ),
        )
    scale = if exact
        one(float_type)
    else
        max(one(float_type), convert(float_type, maximum_magnitude))
    end
    tolerance = converted_absolute + converted_relative * scale
    isfinite(tolerance) || throw(ArgumentError("the combined tolerance is not finite"))
    return tolerance
end

_copclique_exact_scalar(value::Integer) = BigInt(value) // BigInt(1)
function _copclique_exact_scalar(value::Rational)
    return BigInt(numerator(value)) // BigInt(denominator(value))
end
function _copclique_exact_scalar(value::AbstractFloat)
    return rationalize(BigInt, value; tol=zero(value))
end

function _copclique_rational_bits(value::Rational)
    numerator_bits = ndigits(abs(numerator(value)); base=2)
    denominator_bits = ndigits(denominator(value); base=2)
    return max(numerator_bits, denominator_bits)
end

function _copclique_check_bits(value::Rational, limit::Int, label::AbstractString)
    bits = _copclique_rational_bits(value)
    bits <= limit || throw(
        ArgumentError(
            "$label needs $bits exact bits, exceeding max_certificate_bits=$limit"
        ),
    )
    return value
end

function _copclique_unit_simplex(dimension::Int, index::Int)
    vector = fill(BigInt(0) // BigInt(1), dimension)
    vector[index] = BigInt(1) // BigInt(1)
    return vector
end

function _copclique_pair_simplex(dimension::Int, left::Int, right::Int)
    vector = fill(BigInt(0) // BigInt(1), dimension)
    vector[left] = BigInt(1) // BigInt(2)
    vector[right] = BigInt(1) // BigInt(2)
    return vector
end

function _copclique_uniform_simplex(dimension::Int)
    weight = BigInt(1) // BigInt(dimension)
    return fill(weight, dimension)
end

function _copclique_exact_quadratic(
    matrix::AbstractMatrix,
    simplex_vector::AbstractVector{<:Rational},
    max_certificate_bits::Int,
)
    dimension = size(matrix, 1)
    length(simplex_vector) == dimension ||
        throw(DimensionMismatch("simplex witness has the wrong dimension"))
    total = BigInt(0) // BigInt(1)
    for row in 1:dimension
        diagonal = simplex_vector[row]^2 * _copclique_exact_scalar(matrix[row, row])
        total = _copclique_check_bits(
            total + diagonal, max_certificate_bits, "exact quadratic witness"
        )
        for column in (row + 1):dimension
            term =
                2 *
                simplex_vector[row] *
                simplex_vector[column] *
                _copclique_exact_scalar(matrix[row, column])
            total = _copclique_check_bits(
                total + term, max_certificate_bits, "exact quadratic witness"
            )
        end
    end
    return total
end

function _copclique_copositivity_witness(
    simplex_vector::AbstractVector{<:Rational}, value::Rational, source::Symbol
)
    converted = Rational{BigInt}[_copclique_exact_scalar(entry) for entry in simplex_vector]
    read_only = _read_only_plan_array(converted)
    exact_value = _copclique_exact_scalar(value)
    return CopositivityWitness(_VALIDATED_CONSTRUCTOR_TOKEN, read_only, exact_value, source)
end

function _copclique_simplex_from_sphere(point::AbstractVector, max_certificate_bits::Int)
    Base.require_one_based_indexing(point)
    squares = Vector{Rational{BigInt}}(undef, length(point))
    total = BigInt(0) // BigInt(1)
    for index in eachindex(point)
        value = point[index]
        value isa AbstractFloat || throw(
            ArgumentError(
                "hierarchy sample type $(typeof(value)) cannot be rationalized exactly"
            ),
        )
        isfinite(value) ||
            throw(ArgumentError("the hierarchy returned a non-finite sampled point"))
        exact = _copclique_exact_scalar(value)
        square = _copclique_check_bits(
            exact^2, max_certificate_bits, "rationalized sampled coordinate"
        )
        squares[index] = square
        total = _copclique_check_bits(
            total + square, max_certificate_bits, "rationalized sampled norm"
        )
    end
    total > 0 || throw(ArgumentError("the hierarchy returned a zero sampled point"))
    simplex = Vector{Rational{BigInt}}(undef, length(point))
    for index in eachindex(squares)
        simplex[index] = _copclique_check_bits(
            squares[index] / total, max_certificate_bits, "normalized rationalized sample"
        )
    end
    sum(simplex) == 1 || error("internal exact simplex normalization failed")
    return simplex
end

function _copclique_direct_witness(
    matrix::AbstractMatrix, tolerance, max_certificate_bits::Int
)
    dimension = size(matrix, 1)
    exact_tolerance = _copclique_exact_scalar(tolerance)
    boundary = nothing

    for index in 1:dimension
        value = _copclique_check_bits(
            _copclique_exact_scalar(matrix[index, index]),
            max_certificate_bits,
            "coordinate-ray witness",
        )
        if value < -exact_tolerance
            vector = _copclique_unit_simplex(dimension, index)
            return _copclique_copositivity_witness(vector, value, :coordinate_ray), boundary
        elseif value < 0 && boundary === nothing
            boundary = (source=:coordinate_ray, value=value, indices=(index,))
        end
    end

    for left in 1:(dimension - 1), right in (left + 1):dimension
        value = _copclique_check_bits(
            (
                _copclique_exact_scalar(matrix[left, left]) +
                2 * _copclique_exact_scalar(matrix[left, right]) +
                _copclique_exact_scalar(matrix[right, right])
            ) / 4,
            max_certificate_bits,
            "two-coordinate witness",
        )
        if value < -exact_tolerance
            vector = _copclique_pair_simplex(dimension, left, right)
            return _copclique_copositivity_witness(vector, value, :two_coordinate_ray),
            boundary
        elseif value < 0 && boundary === nothing
            boundary = (source=:two_coordinate_ray, value=value, indices=(left, right))
        end
    end

    vector = _copclique_uniform_simplex(dimension)
    value = _copclique_exact_quadratic(matrix, vector, max_certificate_bits)
    if value < -exact_tolerance
        return _copclique_copositivity_witness(vector, value, :uniform_ray), boundary
    elseif value < 0 && boundary === nothing
        boundary = (source=:uniform_ray, value=value, indices=Tuple(1:dimension))
    end
    return nothing, boundary
end

function _copclique_entrywise_nonnegative(matrix::AbstractMatrix)
    minimum_value = matrix[1, 1]
    minimum_index = (1, 1)
    for column in axes(matrix, 2), row in axes(matrix, 1)
        value = matrix[row, column]
        if value < minimum_value
            minimum_value = value
            minimum_index = (row, column)
        end
    end
    return minimum_value >= zero(minimum_value), minimum_value, minimum_index
end

function _copclique_run_hierarchy(
    rng::AbstractRNG,
    polynomial::HomogeneousPolynomial;
    hierarchy::Symbol,
    backend::AbstractOptimizationBackend,
    level::Int,
    sense::Symbol,
    inner_samples::Int,
    allow_densify::Bool,
    limits::OptimizationLimits,
    max_terms,
    max_degree,
    max_polynomial_dimension,
    max_exponent_entries,
    max_dense_entries,
    max_nonzeros,
    max_samples,
    max_full_dimension,
    max_projection_entries,
    max_work,
)
    if hierarchy === :sos
        inner_samples > 0 &&
            !allow_densify &&
            throw(
                ArgumentError(
                    "SOS inner sampling constructs guarded dense hierarchy matrices; " *
                    "pass allow_densify=true explicitly",
                ),
            )
        return polynomial_sos_bounds(
            rng,
            polynomial;
            backend,
            level,
            sense,
            inner_samples,
            limits,
            max_terms,
            max_degree,
            max_dimension=max_polynomial_dimension,
            max_exponent_entries,
            max_dense_entries,
            max_nonzeros,
            max_samples,
            max_full_dimension,
            max_projection_entries,
            max_work,
        )
    end
    return polynomial_bounds(
        rng,
        polynomial;
        level,
        sense,
        inner_samples,
        allow_densify,
        max_terms,
        max_degree,
        max_dimension=max_polynomial_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_samples,
        max_work,
    )
end

function _copclique_hierarchy_bounds(result)
    return result.outer_bound,
    result.inner_bound, result.best_point, result.samples_requested,
    result.samples_evaluated
end

function _copclique_copositivity_status(result, boundary)
    boundary !== nothing && return CopositivityNumericalBoundary
    result isa PolynomialOptimizationResult && return CopositivityHierarchyUnknown
    result.status === OptimizationBackendUnavailable &&
        return CopositivityBackendUnavailable
    result.status === OptimizationLimit && return CopositivityResourceLimit
    result.status in (OptimizationOptimal, OptimizationFeasible) &&
        return CopositivityHierarchyUnknown
    return CopositivityBackendFailure
end

"""
    copositivity_criterion(rng, matrix; hierarchy=:sos,
                          backend=NoOptimizationBackend(), level=0,
                          inner_samples=0, kwargs...)

Analyze whether a finite, real, exactly symmetric matrix `C` is copositive:
`y' * C * y >= 0` for every `y >= 0`.

The routine first checks exact coordinate, pair, and uniform witnesses,
entrywise nonnegativity, and exact PSD sufficiency. A deciding exact branch
returns with `polynomial === nothing`; polynomial and hierarchy limits do not
obstruct that certificate. Otherwise the routine constructs the polynomial
and runs either the solver-neutral SOS hierarchy (`hierarchy=:sos`) or the
guarded dense generalized-eigenvalue hierarchy (`hierarchy=:nosdp`). The RNG
is mandatory, and `inner_samples` is an exact count. No global RNG,
elapsed-time budget, or implicit solver is used.

A numerical lower bound near or above zero is retained but never changed into
`verdict=true`. This corrects the pinned routine's fixed `-1e-9` threshold,
which labels small negative boundaries as copositive. Inputs are never
symmetrized, clipped, normalized, or repaired. Sparse PSD or `:nosdp` work
requires `allow_densify=true`.
"""
function copositivity_criterion(
    rng::AbstractRNG,
    matrix::AbstractMatrix;
    hierarchy=:sos,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    level=0,
    inner_samples=0,
    allow_densify::Bool=false,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    max_matrix_dimension=_COPCLIQUE_DEFAULT_MAX_MATRIX_DIMENSION,
    max_matrix_entries=_COPCLIQUE_DEFAULT_MAX_MATRIX_ENTRIES,
    max_certificate_bits=_COPCLIQUE_DEFAULT_MAX_CERTIFICATE_BITS,
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_polynomial_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_samples=_POLYNOMIAL_DEFAULT_MAX_SAMPLES,
    max_full_dimension=512,
    max_projection_entries=100_000,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    checked_hierarchy = _copclique_hierarchy(hierarchy)
    _copclique_validate_backend(checked_hierarchy, backend)
    checked_level = _copclique_nonnegative_integer(level, "level")
    checked_samples = _copclique_nonnegative_integer(inner_samples, "inner_samples")
    sample_limit = _copclique_optional_limit(max_samples, "max_samples"; allow_zero=true)
    sample_limit !== nothing &&
        checked_samples > sample_limit &&
        throw(
            ArgumentError(
                "inner_samples=$checked_samples exceeds max_samples=$sample_limit"
            ),
        )
    certificate_limit = _copclique_positive_integer(
        max_certificate_bits, "max_certificate_bits"
    )
    _, value_type, exact, maximum_magnitude = _copclique_matrix_contract(
        matrix; graph=false, max_matrix_dimension, max_matrix_entries
    )
    tolerance = _copclique_tolerance(value_type, exact, maximum_magnitude; atol, rtol)

    direct_witness, boundary = _copclique_direct_witness(
        matrix, tolerance, certificate_limit
    )
    if direct_witness !== nothing
        return CopositivityResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            CopositivityCertifiedFalse,
            false,
            true,
            Symbol("exact_", direct_witness.source, "_witness"),
            nothing,
            direct_witness.value,
            :not_computed,
            :attained_exact_witness,
            direct_witness,
            nothing,
            nothing,
            checked_hierarchy,
            checked_level,
            nothing,
            checked_samples,
            0,
            tolerance,
            "an exact rational nonnegative-simplex witness has negative value",
        )
    end

    entrywise, _, _ = _copclique_entrywise_nonnegative(matrix)
    if entrywise
        attained_upper = _copclique_exact_scalar(
            minimum(matrix[index, index] for index in axes(matrix, 1))
        )
        return CopositivityResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            CopositivityCertifiedTrue,
            true,
            true,
            :entrywise_nonnegative,
            BigInt(0) // BigInt(1),
            attained_upper,
            :certified_nonnegative,
            :coordinate_ray,
            nothing,
            nothing,
            nothing,
            checked_hierarchy,
            checked_level,
            nothing,
            checked_samples,
            0,
            tolerance,
            "every matrix entry is nonnegative, which is an exact copositivity certificate",
        )
    end

    psd_diagnostic = if !SparseArrays.issparse(matrix) || allow_densify
        is_positive_semidefinite(matrix; atol=atol, rtol=rtol, allow_densify=allow_densify)
    else
        nothing
    end
    if exact &&
        psd_diagnostic !== nothing &&
        psd_diagnostic.status === MatrixPredicateSatisfied
        attained_upper = minimum(
            _copclique_exact_scalar(matrix[index, index]) for index in axes(matrix, 1)
        )
        return CopositivityResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            CopositivityCertifiedTrue,
            true,
            true,
            :exact_positive_semidefinite,
            BigInt(0) // BigInt(1),
            attained_upper,
            :certified_nonnegative,
            :coordinate_ray,
            nothing,
            psd_diagnostic,
            nothing,
            checked_hierarchy,
            checked_level,
            nothing,
            checked_samples,
            0,
            tolerance,
            "exact positive semidefiniteness is a copositivity certificate",
        )
    end

    polynomial = copositive_polynomial(matrix; max_terms)
    hierarchy_result = _copclique_run_hierarchy(
        rng,
        polynomial;
        hierarchy=checked_hierarchy,
        backend,
        level=checked_level,
        sense=:min,
        inner_samples=checked_samples,
        allow_densify,
        limits,
        max_terms,
        max_degree,
        max_polynomial_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_samples,
        max_full_dimension,
        max_projection_entries,
        max_work,
    )
    lower, upper, best_point, samples_requested, samples_evaluated = _copclique_hierarchy_bounds(
        hierarchy_result
    )

    if best_point !== nothing
        simplex = _copclique_simplex_from_sphere(best_point, certificate_limit)
        exact_value = _copclique_exact_quadratic(matrix, simplex, certificate_limit)
        exact_tolerance = _copclique_exact_scalar(tolerance)
        if exact_value < -exact_tolerance
            sampled_witness = _copclique_copositivity_witness(
                simplex, exact_value, :rationalized_sample
            )
            return CopositivityResult(
                _VALIDATED_CONSTRUCTOR_TOKEN,
                CopositivityCertifiedFalse,
                false,
                true,
                :exact_rationalized_sample_witness,
                lower,
                exact_value,
                hierarchy_result.outer_kind,
                :attained_exact_witness,
                sampled_witness,
                psd_diagnostic,
                polynomial,
                checked_hierarchy,
                checked_level,
                hierarchy_result,
                samples_requested,
                samples_evaluated,
                tolerance,
                "a sampled point was rationalized and independently re-evaluated as an exact negative witness",
            )
        elseif exact_value < 0 && boundary === nothing
            boundary = (source=:rationalized_sample, value=exact_value)
        end
    end

    status = _copclique_copositivity_status(hierarchy_result, boundary)
    message = if status === CopositivityNumericalBoundary
        "negative evidence lies inside the requested tolerance band; copositivity is unknown"
    elseif status === CopositivityBackendUnavailable
        "no SOS backend was supplied; available bounds and diagnostics are retained"
    elseif status === CopositivityResourceLimit
        "the optimizer stopped at a declared limit; no Boolean conclusion is made"
    elseif status === CopositivityBackendFailure
        "the optimizer failed or returned unusable evidence; no Boolean conclusion is made"
    else
        "hierarchy bounds are retained as numerical evidence; they are not an exact copositivity certificate"
    end
    return CopositivityResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        nothing,
        false,
        nothing,
        lower,
        upper,
        hierarchy_result.outer_kind,
        hierarchy_result.inner_kind,
        nothing,
        psd_diagnostic,
        polynomial,
        checked_hierarchy,
        checked_level,
        hierarchy_result,
        samples_requested,
        samples_evaluated,
        tolerance,
        message,
    )
end

function _copclique_graph_integer_copy(adjacency::AbstractMatrix)
    if SparseArrays.issparse(adjacency)
        return sparse(Int.(adjacency))
    end
    return Int.(adjacency)
end

function _copclique_graph_degrees(adjacency::AbstractMatrix)
    dimension = size(adjacency, 1)
    degrees = zeros(Int, dimension)
    for column in 1:dimension
        degree = 0
        for row in 1:dimension
            degree += iszero(adjacency[row, column]) ? 0 : 1
        end
        degrees[column] = degree
    end
    return degrees
end

function _copclique_check_greedy_work(dimension::Int, max_greedy_work)
    work = BigInt(dimension)^3 + 5 * BigInt(dimension)^2
    work_limit = _copclique_optional_limit(
        max_greedy_work, "max_greedy_work"; allow_zero=true
    )
    work_limit !== nothing &&
        work > work_limit &&
        throw(
            ArgumentError(
                "deterministic graph certificates need at most $work adjacency checks, " *
                "exceeding max_greedy_work=$work_limit",
            ),
        )
    return work
end

function _copclique_best_greedy_clique(adjacency::AbstractMatrix, degrees::Vector{Int})
    dimension = size(adjacency, 1)
    order = sort(collect(1:dimension); by=vertex -> (-degrees[vertex], vertex))
    best = Int[]
    for start in 1:dimension
        candidate = Int[]
        for offset in 0:(dimension - 1)
            vertex = order[mod1(start + offset, dimension)]
            all(other -> !iszero(adjacency[vertex, other]), candidate) &&
                push!(candidate, vertex)
        end
        if length(candidate) > length(best)
            best = candidate
        end
    end
    sort!(best)
    return best, order
end

function _copclique_greedy_coloring_upper(adjacency::AbstractMatrix, order::Vector{Int})
    dimension = size(adjacency, 1)
    colors = zeros(Int, dimension)
    forbidden = falses(dimension)
    color_count = 0
    for vertex in order
        fill!(forbidden, false)
        for neighbor in 1:dimension
            if !iszero(adjacency[vertex, neighbor]) && colors[neighbor] > 0
                forbidden[colors[neighbor]] = true
            end
        end
        color = findfirst(!, forbidden)
        color === nothing && error("internal greedy-color selection failed")
        colors[vertex] = color
        color_count = max(color_count, color)
    end
    return color_count, colors
end

function _copclique_edge_upper(edge_count::BigInt, dimension::Int)
    root = isqrt(8 * edge_count + 1)
    return min(dimension, Int((root + 1) ÷ 2))
end

function _copclique_upper_certificate(
    adjacency::AbstractMatrix,
    degrees::Vector{Int},
    edge_count::BigInt,
    best_clique::Vector{Int},
    colors::Vector{Int},
    edge_upper::Int,
    degree_upper::Int,
    coloring_upper::Int,
    structural_upper::Int,
)
    dimension = size(adjacency, 1)
    length(degrees) == dimension ||
        error("internal graph degree certificate has the wrong dimension")
    length(colors) == dimension ||
        error("internal graph coloring certificate has the wrong dimension")
    issorted(best_clique) && allunique(best_clique) ||
        error("internal greedy clique certificate is not sorted and unique")
    all(vertex -> 1 <= vertex <= dimension, best_clique) ||
        error("internal greedy clique certificate has an invalid vertex")
    for position in eachindex(best_clique)
        left = best_clique[position]
        for right in best_clique[(position + 1):end]
            !iszero(adjacency[left, right]) ||
                error("internal greedy clique certificate contains a non-edge")
        end
    end

    checked_degree_sum = BigInt(0)
    for vertex in 1:dimension
        checked_degree = 0
        for neighbor in 1:dimension
            checked_degree += iszero(adjacency[vertex, neighbor]) ? 0 : 1
        end
        checked_degree == degrees[vertex] ||
            error("internal maximum-degree certificate is inconsistent")
        checked_degree_sum += checked_degree
    end
    checked_degree_sum == 2 * edge_count ||
        error("internal edge-count certificate is inconsistent")

    all(color -> color >= 1, colors) ||
        error("internal graph coloring uses a nonpositive color")
    maximum(colors) == coloring_upper ||
        error("internal graph coloring does not attain its recorded bound")
    for left in 1:(dimension - 1), right in (left + 1):dimension
        if !iszero(adjacency[left, right])
            colors[left] != colors[right] ||
                error("internal graph coloring assigns one color to adjacent vertices")
        end
    end

    maximum_degree = maximum(degrees)
    edge_upper == _copclique_edge_upper(edge_count, dimension) ||
        error("internal edge-count upper bound is inconsistent")
    degree_upper == min(dimension, maximum_degree + 1) ||
        error("internal maximum-degree upper bound is inconsistent")
    structural_upper == min(edge_upper, degree_upper, coloring_upper) ||
        error("internal combined graph upper bound is inconsistent")
    attaining_kinds = Symbol[]
    edge_upper == structural_upper && push!(attaining_kinds, :edge_count)
    degree_upper == structural_upper && push!(attaining_kinds, :maximum_degree)
    coloring_upper == structural_upper && push!(attaining_kinds, :greedy_coloring)
    isempty(attaining_kinds) && error("internal graph upper bound has no certificate")
    return (
        edge_count_bound=edge_upper,
        maximum_degree=maximum_degree,
        maximum_degree_bound=degree_upper,
        greedy_coloring_bound=coloring_upper,
        coloring=_read_only_plan_array(colors),
        attaining_kinds=Tuple(attaining_kinds),
    )
end

function _copclique_numerical_upper_candidate(value, tolerance, dimension::Int)
    value === nothing && return nothing
    padded = value + tolerance
    padded >= one(padded) && return dimension
    padded <= zero(padded) && return 1
    reciprocal = inv(one(padded) - padded)
    isfinite(reciprocal) || return dimension
    reciprocal >= dimension && return dimension
    return clamp(floor(Int, reciprocal), 1, dimension)
end

function _copclique_motzkin_straus_witness(
    adjacency::AbstractMatrix, point, max_certificate_bits::Int
)
    point === nothing && return nothing
    simplex = _copclique_simplex_from_sphere(point, max_certificate_bits)
    value = _copclique_exact_quadratic(adjacency, simplex, max_certificate_bits)
    value < 1 ||
        throw(ArgumentError("the rationalized Motzkin--Straus value is not below one"))
    reciprocal = inv(1 - value)
    implied = cld(numerator(reciprocal), denominator(reciprocal))
    dimension = size(adjacency, 1)
    1 <= implied <= dimension || throw(
        ArgumentError(
            "the rationalized Motzkin--Straus lower bound $implied is outside 1:$dimension",
        ),
    )
    read_only = _read_only_plan_array(Rational{BigInt}.(simplex))
    return MotzkinStrausWitness(
        _VALIDATED_CONSTRUCTOR_TOKEN, read_only, value, Int(implied)
    )
end

function _copclique_clique_status(hierarchy_result)
    hierarchy_result isa PolynomialOptimizationResult &&
        return CliqueNumberNumericalHierarchy
    hierarchy_result.status === OptimizationBackendUnavailable &&
        return CliqueNumberBackendUnavailable
    hierarchy_result.status === OptimizationLimit && return CliqueNumberResourceLimit
    hierarchy_result.status in (OptimizationOptimal, OptimizationFeasible) &&
        return CliqueNumberNumericalHierarchy
    return CliqueNumberBackendFailure
end

"""
    clique_number_bounds(rng, adjacency; hierarchy=:sos,
                         backend=NoOptimizationBackend(), level=0,
                         inner_samples=0, kwargs...)

Return certified integer lower and upper bounds on the clique number of a
nonempty simple undirected graph. The adjacency matrix must be finite, exactly
symmetric, zero on the diagonal, and contain only exact zero/one values.

The lower certificate is a deterministic greedy clique, optionally improved
by exact rational re-evaluation of an attained Motzkin--Straus point. Exact
edge-count, degree, and greedy-coloring arguments provide the upper bound.
The result retains the proper coloring and every component upper bound so the
certificate can be checked against the input graph. If the graph certificates
already meet, the routine returns with `polynomial === nothing`. Otherwise
the selected polynomial hierarchy is retained separately; a floating outer
value only produces `uncertified_upper_candidate`.

This separation corrects two unsafe pinned behaviors: elapsed-time-dependent
global-RNG sampling and fixed decimal offsets before integer rounding.
The edgeless graph is handled exactly as clique number one. Inputs are never
rounded, symmetrized, binarized, or otherwise repaired.
"""
function clique_number_bounds(
    rng::AbstractRNG,
    adjacency::AbstractMatrix;
    hierarchy=:sos,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    level=0,
    inner_samples=0,
    allow_densify::Bool=false,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    max_matrix_dimension=_COPCLIQUE_DEFAULT_MAX_MATRIX_DIMENSION,
    max_matrix_entries=_COPCLIQUE_DEFAULT_MAX_MATRIX_ENTRIES,
    max_greedy_work=_COPCLIQUE_DEFAULT_MAX_GREEDY_WORK,
    max_certificate_bits=_COPCLIQUE_DEFAULT_MAX_CERTIFICATE_BITS,
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_polynomial_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_samples=_POLYNOMIAL_DEFAULT_MAX_SAMPLES,
    max_full_dimension=512,
    max_projection_entries=100_000,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    checked_hierarchy = _copclique_hierarchy(hierarchy)
    _copclique_validate_backend(checked_hierarchy, backend)
    checked_level = _copclique_nonnegative_integer(level, "level")
    checked_samples = _copclique_nonnegative_integer(inner_samples, "inner_samples")
    sample_limit = _copclique_optional_limit(max_samples, "max_samples"; allow_zero=true)
    sample_limit !== nothing &&
        checked_samples > sample_limit &&
        throw(
            ArgumentError(
                "inner_samples=$checked_samples exceeds max_samples=$sample_limit"
            ),
        )
    certificate_limit = _copclique_positive_integer(
        max_certificate_bits, "max_certificate_bits"
    )
    dimension, value_type, exact, _ = _copclique_matrix_contract(
        adjacency; graph=true, max_matrix_dimension, max_matrix_entries
    )
    exact ||
        value_type <: AbstractFloat ||
        throw(ArgumentError("unsupported adjacency element type $value_type"))
    tolerance = _copclique_tolerance(value_type, exact, one(value_type); atol, rtol)
    _copclique_check_greedy_work(dimension, max_greedy_work)
    integer_adjacency = _copclique_graph_integer_copy(adjacency)
    degrees = _copclique_graph_degrees(integer_adjacency)
    degree_sum = sum(BigInt, degrees)
    iseven(degree_sum) || error("internal graph degree sum is odd")
    edge_count = degree_sum ÷ 2
    best_clique, order = _copclique_best_greedy_clique(integer_adjacency, degrees)
    coloring_upper, colors = _copclique_greedy_coloring_upper(integer_adjacency, order)
    edge_upper = _copclique_edge_upper(edge_count, dimension)
    degree_upper = min(dimension, maximum(degrees) + 1)
    structural_upper = min(edge_upper, degree_upper, coloring_upper)
    structural_lower = length(best_clique)
    structural_lower <= structural_upper ||
        error("internal graph certificates produced inconsistent bounds")
    upper_certificate = _copclique_upper_certificate(
        integer_adjacency,
        degrees,
        edge_count,
        best_clique,
        colors,
        edge_upper,
        degree_upper,
        coloring_upper,
        structural_upper,
    )
    upper_kinds = upper_certificate.attaining_kinds
    read_only_best_clique = _read_only_plan_array(best_clique)

    if structural_lower == structural_upper
        return CliqueNumberResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            CliqueNumberExact,
            structural_lower,
            structural_upper,
            true,
            true,
            :greedy_clique,
            upper_kinds,
            upper_certificate,
            read_only_best_clique,
            dimension,
            edge_count,
            nothing,
            checked_hierarchy,
            checked_level,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            checked_samples,
            0,
            tolerance,
            "deterministic graph certificates determine the clique number exactly; no RNG was consumed",
        )
    end

    polynomial = copositive_polynomial(integer_adjacency; max_terms)
    hierarchy_result = _copclique_run_hierarchy(
        rng,
        polynomial;
        hierarchy=checked_hierarchy,
        backend,
        level=checked_level,
        sense=:max,
        inner_samples=checked_samples,
        allow_densify,
        limits,
        max_terms,
        max_degree,
        max_polynomial_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_samples,
        max_full_dimension,
        max_projection_entries,
        max_work,
    )
    continuous_upper, continuous_lower, best_point, samples_requested, samples_evaluated = _copclique_hierarchy_bounds(
        hierarchy_result
    )
    numerical_upper = _copclique_numerical_upper_candidate(
        continuous_upper, tolerance, dimension
    )
    motzkin_witness = _copclique_motzkin_straus_witness(
        integer_adjacency, best_point, certificate_limit
    )
    sampled_lower =
        motzkin_witness === nothing ? structural_lower : motzkin_witness.implied_lower_bound
    certified_lower = max(structural_lower, sampled_lower)
    lower_kind =
        certified_lower > structural_lower ? :motzkin_straus_witness : :greedy_clique
    certified_lower <= structural_upper ||
        error("independent clique certificates produced inconsistent bounds")
    exact_result = certified_lower == structural_upper
    status = exact_result ? CliqueNumberExact : _copclique_clique_status(hierarchy_result)
    message = if exact_result
        "exact graph certificates and an attained Motzkin--Straus witness determine the clique number"
    elseif status === CliqueNumberBackendUnavailable
        "no SOS backend was supplied; the returned discrete graph bounds remain certified"
    elseif status === CliqueNumberResourceLimit
        "the optimizer stopped at a declared limit; the returned discrete graph bounds remain certified"
    elseif status === CliqueNumberBackendFailure
        "the optimizer failed or returned unusable evidence; the returned discrete graph bounds remain certified"
    else
        "the discrete interval is certified; the floating hierarchy upper candidate is retained but not promoted to a certificate"
    end
    return CliqueNumberResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        certified_lower,
        structural_upper,
        exact_result,
        true,
        lower_kind,
        upper_kinds,
        upper_certificate,
        read_only_best_clique,
        dimension,
        edge_count,
        polynomial,
        checked_hierarchy,
        checked_level,
        hierarchy_result,
        continuous_lower,
        continuous_upper,
        numerical_upper,
        motzkin_witness,
        samples_requested,
        samples_evaluated,
        tolerance,
        message,
    )
end

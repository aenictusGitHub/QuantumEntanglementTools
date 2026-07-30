# Source-informed independent Julia implementation based on the specification
# and QETLAB FilterNormalForm.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export FilterNormalFormResult, filter_normal_form

const _FILTER_NORMAL_FORM_DEFAULT_MAX_ENTRIES = 10_000_000
const _FILTER_NORMAL_FORM_DEFAULT_MAX_WORK = 1_000_000_000

"""
    FilterNormalFormResult

Structured result returned by [`filter_normal_form`](@ref).

For `status == :converged`, let `τ = normal_form_trace`,
`n = input_dimension`, and `K = tensor_product(local_filters...)`. Then

```math
\\mathtt{filtered\\_operator} \\approx K\\rho K^\\dagger
\\approx \\frac{\\tau I +
  \\sum_k \\mathtt{coefficients}_k
  \\mathtt{left\\_operators}_k \\otimes
  \\mathtt{right\\_operators}_k}{n}.
```

The local operators are Frobenius-orthonormal and Hermitian. All coefficients
from the thin operator-Schmidt decomposition are retained, including
numerical zeros. `coefficient_numerical_rank` and `coefficient_threshold`
report a numerical-rank diagnostic without silently deleting terms.

`status` is one of `:converged`, `:max_iterations`, `:singular_marginal`,
`:ill_conditioned`, `:work_limit`, or `:numerical_failure`. The first five
mirror the bounded Sinkhorn stage. `:numerical_failure` also covers a derived
Hermiticity or reconstruction check that crossed its explicit tolerance.
Only `:converged` sets `converged=true`.

`input_numerical_rank`, its threshold, extremal input eigenvalues, and
`input_numerically_full_rank` make rank deficiency explicit. Rank deficiency
is diagnostic rather than an automatic failure: some low-rank operators
already have balanced nonsingular marginals. The nested `sinkhorn_result`,
filter condition data, residual history, and failed subsystem remain
available through the Sinkhorn result whenever that stage started.

No input or returned operator is normalized, symmetrized, clipped, or
otherwise repaired. A derived centered work matrix is averaged with its
adjoint only after its Hermiticity residual is checked against and recorded
with `decomposition_hermiticity_tolerance`.
"""
struct FilterNormalFormResult{
    FO,CO,C,LO,RO,FS,LF,RF,S,FT,RT,EMIN,EMAX,CT,DHR,DHT,NRR,NRT,FIR,FIT
}
    filtered_operator::FO
    centered_operator::CO
    coefficients::C
    left_operators::LO
    right_operators::RO
    local_filters::FS
    left_filter::LF
    right_filter::RF
    sinkhorn_result::S
    status::Symbol
    converged::Bool
    input_dimension::Int
    normal_form_trace::FT
    input_numerical_rank::Union{Nothing,Int}
    input_rank_threshold::RT
    input_minimum_eigenvalue::EMIN
    input_maximum_eigenvalue::EMAX
    input_numerically_full_rank::Union{Nothing,Bool}
    coefficient_numerical_rank::Union{Nothing,Int}
    coefficient_threshold::CT
    decomposition_hermiticity_residual::DHR
    decomposition_hermiticity_tolerance::DHT
    normal_form_reconstruction_residual::NRR
    normal_form_reconstruction_tolerance::NRT
    filter_identity_residual::FIR
    filter_identity_tolerance::FIT
    work_used::BigInt
    work_limit::Union{Nothing,BigInt}
    estimated_postprocessing_work::BigInt
    required_next_stage_work::BigInt
    failed_subsystem::Union{Nothing,Int}
    message::String
end

function Base.show(io::IO, result::FilterNormalFormResult)
    return print(
        io,
        "FilterNormalFormResult(status=",
        result.status,
        ", input_rank=",
        result.input_numerical_rank,
        "/",
        result.input_dimension,
        ", coefficient_rank=",
        result.coefficient_numerical_rank,
        ", work_used=",
        result.work_used,
        ")",
    )
end

function _filter_normal_form_result(;
    filtered_operator=nothing,
    centered_operator=nothing,
    coefficients=nothing,
    left_operators=nothing,
    right_operators=nothing,
    local_filters=nothing,
    left_filter=nothing,
    right_filter=nothing,
    sinkhorn_result=nothing,
    status,
    input_dimension,
    normal_form_trace=nothing,
    input_numerical_rank=nothing,
    input_rank_threshold=nothing,
    input_minimum_eigenvalue=nothing,
    input_maximum_eigenvalue=nothing,
    input_numerically_full_rank=nothing,
    coefficient_numerical_rank=nothing,
    coefficient_threshold=nothing,
    decomposition_hermiticity_residual=nothing,
    decomposition_hermiticity_tolerance=nothing,
    normal_form_reconstruction_residual=nothing,
    normal_form_reconstruction_tolerance=nothing,
    filter_identity_residual=nothing,
    filter_identity_tolerance=nothing,
    work_used=BigInt(0),
    work_limit=nothing,
    estimated_postprocessing_work=BigInt(0),
    required_next_stage_work=BigInt(0),
    failed_subsystem=nothing,
    message,
)
    return FilterNormalFormResult(
        filtered_operator,
        centered_operator,
        coefficients,
        left_operators,
        right_operators,
        local_filters,
        left_filter,
        right_filter,
        sinkhorn_result,
        status,
        status === :converged,
        input_dimension,
        normal_form_trace,
        input_numerical_rank,
        input_rank_threshold,
        input_minimum_eigenvalue,
        input_maximum_eigenvalue,
        input_numerically_full_rank,
        coefficient_numerical_rank,
        coefficient_threshold,
        decomposition_hermiticity_residual,
        decomposition_hermiticity_tolerance,
        normal_form_reconstruction_residual,
        normal_form_reconstruction_tolerance,
        filter_identity_residual,
        filter_identity_tolerance,
        BigInt(work_used),
        work_limit,
        BigInt(estimated_postprocessing_work),
        BigInt(required_next_stage_work),
        failed_subsystem,
        message,
    )
end

function _filter_normal_form_tolerance(value, name::AbstractString, ::Type{R}) where {R}
    checked = _tierd_validate_tolerance(value, name)
    checked === nothing && return nothing
    converted = try
        R(checked)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$(repr(value)) cannot be represented as $R"))
    end
    isfinite(converted) ||
        throw(ArgumentError("$name=$(repr(value)) is not finite when represented as $R"))
    return converted
end

function _filter_normal_form_tolerance_pair(
    ::Type{R}, atol, rtol, atol_name, rtol_name, default_rtol
) where {R}
    absolute = _filter_normal_form_tolerance(atol, atol_name, R)
    relative = _filter_normal_form_tolerance(rtol, rtol_name, R)
    return (
        absolute === nothing ? zero(R) : absolute,
        relative === nothing ? R(default_rtol) : relative,
    )
end

function _filter_normal_form_budget(
    dimension::Int, layout::SubsystemLayout, max_entries, max_work
)
    entry_limit = _operator_sinkhorn_limit(max_entries, "max_entries")
    work_limit = _operator_sinkhorn_limit(max_work, "max_work")
    left_operator_dimension = BigInt(layout[1])^2
    right_operator_dimension = BigInt(layout[2])^2
    term_count = min(left_operator_dimension, right_operator_dimension)
    full_entries = BigInt(dimension)^2

    # This guards the dense input/Sinkhorn/centered matrices, both explicit
    # Hermitian operator bases used by the decomposition, and its coordinate
    # matrix. It is deliberately conservative and documented as a total
    # materialization estimate rather than a peak-memory claim.
    estimated_dense_entries =
        4 * full_entries + left_operator_dimension^2 + right_operator_dimension^2
    entry_limit !== nothing &&
        estimated_dense_entries > entry_limit &&
        throw(
            ArgumentError(
                "filter_normal_form has a conservative dense materialization " *
                "estimate of $estimated_dense_entries entries, exceeding " *
                "max_entries=$entry_limit",
            ),
        )

    sinkhorn_budget = _operator_sinkhorn_budget(dimension, layout, nothing, nothing)
    rank_work = BigInt(dimension)^3
    svd_work = full_entries * term_count
    reconstruction_work = full_entries * term_count
    filter_identity_work = 2 * BigInt(dimension)^3
    postprocessing_work =
        2 * full_entries + svd_work + reconstruction_work + filter_identity_work
    return (;
        work_limit, sinkhorn_budget, rank_work, postprocessing_work, estimated_dense_entries
    )
end

function _filter_normal_form_sinkhorn_failure(
    sinkhorn,
    dimension,
    rank_data,
    rank_threshold,
    work_used,
    work_limit,
    postprocessing_work,
    required_next_stage_work,
)
    return _filter_normal_form_result(;
        filtered_operator=sinkhorn.scaled_operator,
        local_filters=sinkhorn.local_filters,
        left_filter=sinkhorn.left_filter,
        right_filter=sinkhorn.right_filter,
        sinkhorn_result=sinkhorn,
        status=sinkhorn.status,
        input_dimension=dimension,
        normal_form_trace=tr(sinkhorn.scaled_operator),
        input_numerical_rank=rank_data.numerical_rank,
        input_rank_threshold=rank_threshold,
        input_minimum_eigenvalue=rank_data.minimum_eigenvalue,
        input_maximum_eigenvalue=rank_data.maximum_eigenvalue,
        input_numerically_full_rank=rank_data.numerical_rank == dimension,
        work_used=work_used,
        work_limit=work_limit,
        estimated_postprocessing_work=postprocessing_work,
        required_next_stage_work=required_next_stage_work,
        failed_subsystem=sinkhorn.failed_subsystem,
        message="operator Sinkhorn stage ended with status $(sinkhorn.status): " *
                sinkhorn.message,
    )
end

"""
    filter_normal_form(
        rho, dims;
        balance_atol=0,
        balance_rtol=nothing,
        rank_atol=0,
        rank_rtol=nothing,
        coefficient_atol=0,
        coefficient_rtol=nothing,
        verification_atol=0,
        verification_rtol=nothing,
        max_iterations=1000,
        max_condition_number=nothing,
        allow_densify=false,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    ) -> FilterNormalFormResult

Compute a bounded bipartite filter normal form. `rho` must be a finite,
exactly Hermitian positive-semidefinite matrix with positive trace and
`prod(dims) == size(rho, 1)`. Exactly two subsystems are required.

The balancing stage delegates to [`operator_sinkhorn`](@ref). Its aggregate
marginal threshold is
`balance_atol + balance_rtol * sum(1 / sqrt(dims[j]))`, with
`balance_rtol = sqrt(eps(R))` by default. Iteration, condition, entry, and
conservative work bounds are mandatory by default and are never inferred
from warnings or global state.

The input numerical-rank threshold is
`rank_atol + rank_rtol * maximum(eigvals(rho))`, where `rank_rtol` defaults
to `size(rho, 1) * eps(R)`. The coefficient-rank threshold has the same form
and defaults to `max(dims .^ 2) * eps(R)`, matching the scale of QETLAB's
default operator-Schmidt truncation while retaining all terms in the native
result.

The normal-form and filter identities are independently checked with
`verification_atol + verification_rtol * max(1, scale)`;
`verification_rtol` defaults to `128 * size(rho, 1) * eps(R)`. Crossing
either bound returns `:numerical_failure` with the measured residual rather
than returning an unchecked decomposition.

Sparse input requires `allow_densify=true`; output is dense. The core accepts
only `Float32`, `Float64`, `ComplexF32`, and `ComplexF64`, preserving that
element type. No matrix inverse is formed and no precision-changing
conversion, input normalization, spectral clipping, or output repair occurs.
"""
function filter_normal_form(
    rho::AbstractMatrix{<:Number},
    dims;
    balance_atol=0,
    balance_rtol=nothing,
    rank_atol=0,
    rank_rtol=nothing,
    coefficient_atol=0,
    coefficient_rtol=nothing,
    verification_atol=0,
    verification_rtol=nothing,
    max_iterations=_OPERATOR_SINKHORN_DEFAULT_MAX_ITERATIONS,
    max_condition_number=nothing,
    allow_densify::Bool=false,
    max_entries=_FILTER_NORMAL_FORM_DEFAULT_MAX_ENTRIES,
    max_work=_FILTER_NORMAL_FORM_DEFAULT_MAX_WORK,
)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("filter_normal_form requires a square matrix"))
    size(rho, 1) > 0 || throw(ArgumentError("rho must have positive dimension"))
    layout = _as_layout(dims)
    length(layout) == 2 ||
        throw(ArgumentError("filter_normal_form requires exactly two subsystems"))
    _validate_matrix_dimension(rho, layout)
    checked_iterations = _operator_sinkhorn_nonnegative_integer(
        max_iterations, "max_iterations"
    )
    max_condition_number === nothing ||
        _operator_sinkhorn_positive_real(max_condition_number, "max_condition_number")

    budget = _filter_normal_form_budget(size(rho, 1), layout, max_entries, max_work)
    dense, real_type = _tierd_dense_matrix(
        rho; allow_densify=allow_densify, operation="filter_normal_form"
    )
    ishermitian(dense) || throw(
        ArgumentError(
            "filter_normal_form requires an exactly Hermitian rho; the input " *
            "is never symmetrized",
        ),
    )
    input_trace = real(tr(dense))
    isfinite(input_trace) ||
        throw(ArgumentError("rho has a non-finite trace after accumulation"))
    input_trace > zero(input_trace) ||
        throw(DomainError(input_trace, "rho must have positive trace"))

    checked_balance_atol, checked_balance_rtol = _filter_normal_form_tolerance_pair(
        real_type,
        balance_atol,
        balance_rtol,
        "balance_atol",
        "balance_rtol",
        sqrt(eps(real_type)),
    )
    checked_rank_atol, checked_rank_rtol = _filter_normal_form_tolerance_pair(
        real_type,
        rank_atol,
        rank_rtol,
        "rank_atol",
        "rank_rtol",
        real_type(size(rho, 1)) * eps(real_type),
    )
    checked_coefficient_atol, checked_coefficient_rtol = _filter_normal_form_tolerance_pair(
        real_type,
        coefficient_atol,
        coefficient_rtol,
        "coefficient_atol",
        "coefficient_rtol",
        real_type(max(layout[1]^2, layout[2]^2)) * eps(real_type),
    )
    checked_verification_atol, checked_verification_rtol = _filter_normal_form_tolerance_pair(
        real_type,
        verification_atol,
        verification_rtol,
        "verification_atol",
        "verification_rtol",
        128 * real_type(size(rho, 1)) * eps(real_type),
    )

    minimum_start_work = budget.rank_work + budget.sinkhorn_budget.setup_work
    if budget.work_limit !== nothing && minimum_start_work > budget.work_limit
        return _filter_normal_form_result(;
            status=:work_limit,
            input_dimension=size(rho, 1),
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            required_next_stage_work=minimum_start_work,
            message="input rank validation plus the initial Sinkhorn snapshot " *
                    "requires estimated work $minimum_start_work, exceeding " *
                    "max_work=$(budget.work_limit)",
        )
    end

    input_decomposition = try
        eigen(Hermitian(dense))
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        return _filter_normal_form_result(;
            status=:numerical_failure,
            input_dimension=size(rho, 1),
            work_used=budget.rank_work,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            message="the checked Hermitian input eigendecomposition failed",
        )
    end
    minimum_input_eigenvalue = minimum(input_decomposition.values)
    maximum_input_eigenvalue = maximum(input_decomposition.values)
    if !all(isfinite, input_decomposition.values)
        return _filter_normal_form_result(;
            status=:numerical_failure,
            input_dimension=size(rho, 1),
            input_minimum_eigenvalue=minimum_input_eigenvalue,
            input_maximum_eigenvalue=maximum_input_eigenvalue,
            work_used=budget.rank_work,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            message="the checked Hermitian input eigendecomposition returned a " *
                    "non-finite eigenvalue",
        )
    end
    minimum_input_eigenvalue >= zero(real_type) || throw(
        DomainError(
            minimum_input_eigenvalue,
            "rho must be positive semidefinite; negative eigenvalues are never clipped",
        ),
    )
    input_rank_threshold = checked_rank_atol + checked_rank_rtol * maximum_input_eigenvalue
    input_numerical_rank = count(>(input_rank_threshold), input_decomposition.values)
    rank_data = (
        numerical_rank=input_numerical_rank,
        minimum_eigenvalue=minimum_input_eigenvalue,
        maximum_eigenvalue=maximum_input_eigenvalue,
    )

    sinkhorn_work_limit = if budget.work_limit === nothing
        nothing
    else
        budget.work_limit - budget.rank_work
    end
    sinkhorn = operator_sinkhorn(
        dense,
        layout;
        atol=checked_balance_atol,
        rtol=checked_balance_rtol,
        max_iterations=checked_iterations,
        max_condition_number=max_condition_number,
        allow_densify=false,
        max_entries=max_entries,
        max_work=sinkhorn_work_limit,
    )
    work_used = budget.rank_work + sinkhorn.work_used
    if !sinkhorn.converged
        required_next_stage_work =
            sinkhorn.status === :work_limit ? budget.sinkhorn_budget.sweep_work : BigInt(0)
        return _filter_normal_form_sinkhorn_failure(
            sinkhorn,
            size(rho, 1),
            rank_data,
            input_rank_threshold,
            work_used,
            budget.work_limit,
            budget.postprocessing_work,
            required_next_stage_work,
        )
    end

    if budget.work_limit !== nothing &&
        work_used + budget.postprocessing_work > budget.work_limit
        return _filter_normal_form_result(;
            filtered_operator=sinkhorn.scaled_operator,
            local_filters=sinkhorn.local_filters,
            left_filter=sinkhorn.left_filter,
            right_filter=sinkhorn.right_filter,
            sinkhorn_result=sinkhorn,
            status=:work_limit,
            input_dimension=size(rho, 1),
            normal_form_trace=tr(sinkhorn.scaled_operator),
            input_numerical_rank=input_numerical_rank,
            input_rank_threshold=input_rank_threshold,
            input_minimum_eigenvalue=minimum_input_eigenvalue,
            input_maximum_eigenvalue=maximum_input_eigenvalue,
            input_numerically_full_rank=input_numerical_rank == size(rho, 1),
            work_used=work_used,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            required_next_stage_work=budget.postprocessing_work,
            message="normal-form decomposition and identity verification require " *
                    "estimated work $(budget.postprocessing_work), exceeding the " *
                    "remaining max_work budget",
        )
    end

    filtered_operator = sinkhorn.scaled_operator
    dimension = size(filtered_operator, 1)
    normal_form_trace = real(tr(filtered_operator))
    identity_matrix = Matrix{eltype(filtered_operator)}(I, dimension, dimension)
    centered_operator =
        filtered_operator - (normal_form_trace / real_type(dimension)) * identity_matrix
    decomposition_hermiticity_residual = maximum(
        abs, centered_operator - adjoint(centered_operator); init=zero(real_type)
    )
    centered_scale = maximum(abs, centered_operator; init=zero(real_type))
    decomposition_hermiticity_tolerance =
        checked_verification_atol +
        checked_verification_rtol * max(one(real_type), centered_scale)
    if !isfinite(normal_form_trace) ||
        !isfinite(decomposition_hermiticity_residual) ||
        !isfinite(decomposition_hermiticity_tolerance) ||
        decomposition_hermiticity_residual > decomposition_hermiticity_tolerance
        return _filter_normal_form_result(;
            filtered_operator=filtered_operator,
            centered_operator=centered_operator,
            local_filters=sinkhorn.local_filters,
            left_filter=sinkhorn.left_filter,
            right_filter=sinkhorn.right_filter,
            sinkhorn_result=sinkhorn,
            status=:numerical_failure,
            input_dimension=dimension,
            normal_form_trace=normal_form_trace,
            input_numerical_rank=input_numerical_rank,
            input_rank_threshold=input_rank_threshold,
            input_minimum_eigenvalue=minimum_input_eigenvalue,
            input_maximum_eigenvalue=maximum_input_eigenvalue,
            input_numerically_full_rank=input_numerical_rank == dimension,
            decomposition_hermiticity_residual=decomposition_hermiticity_residual,
            decomposition_hermiticity_tolerance=decomposition_hermiticity_tolerance,
            work_used=work_used,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            message="the derived centered operator has Hermiticity residual " *
                    "$decomposition_hermiticity_residual above the checked " *
                    "work-matrix tolerance $decomposition_hermiticity_tolerance",
        )
    end

    # This projection is confined to a derived decomposition work matrix and
    # its measured residual is returned above. The input, filtered operator,
    # and centered operator remain unmodified.
    decomposition_work_matrix =
        (centered_operator + adjoint(centered_operator)) / real_type(2)
    decomposition = try
        operator_schmidt_decomposition(
            decomposition_work_matrix,
            layout;
            allow_densify=false,
            hermitian_factors=true,
        )
    catch err
        coordinate_failure =
            err isa ErrorException &&
            startswith(err.msg, "Hermitian-basis coordinates have imaginary residual")
        (err isa LinearAlgebra.LAPACKException || coordinate_failure) || rethrow()
        return _filter_normal_form_result(;
            filtered_operator=filtered_operator,
            centered_operator=centered_operator,
            local_filters=sinkhorn.local_filters,
            left_filter=sinkhorn.left_filter,
            right_filter=sinkhorn.right_filter,
            sinkhorn_result=sinkhorn,
            status=:numerical_failure,
            input_dimension=dimension,
            normal_form_trace=normal_form_trace,
            input_numerical_rank=input_numerical_rank,
            input_rank_threshold=input_rank_threshold,
            input_minimum_eigenvalue=minimum_input_eigenvalue,
            input_maximum_eigenvalue=maximum_input_eigenvalue,
            input_numerically_full_rank=input_numerical_rank == dimension,
            decomposition_hermiticity_residual=decomposition_hermiticity_residual,
            decomposition_hermiticity_tolerance=decomposition_hermiticity_tolerance,
            work_used=work_used + budget.postprocessing_work,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            message="the checked Hermitian operator-Schmidt decomposition failed",
        )
    end

    coefficients = real_type(dimension) .* decomposition.coefficients
    coefficient_scale = maximum(coefficients; init=zero(real_type))
    coefficient_threshold =
        checked_coefficient_atol + checked_coefficient_rtol * coefficient_scale
    coefficient_numerical_rank = count(>(coefficient_threshold), coefficients)
    tensor_component = tensor_sum(
        decomposition.left_factors, decomposition.right_factors; weights=coefficients
    )
    normal_form_reconstruction =
        (normal_form_trace * identity_matrix + tensor_component) / real_type(dimension)
    normal_form_reconstruction_residual = norm(
        filtered_operator - normal_form_reconstruction
    )
    normal_form_reconstruction_tolerance =
        checked_verification_atol +
        checked_verification_rtol * max(one(real_type), norm(filtered_operator))

    full_filter = tensor_product(sinkhorn.local_filters...)
    filter_identity_operator = full_filter * dense * adjoint(full_filter)
    filter_identity_residual = norm(filtered_operator - filter_identity_operator)
    filter_identity_tolerance =
        checked_verification_atol +
        checked_verification_rtol *
        max(one(real_type), norm(filtered_operator), norm(filter_identity_operator))
    work_used += budget.postprocessing_work

    decomposition_is_finite =
        all(isfinite, coefficients) &&
        all(factor -> all(isfinite, factor), decomposition.left_factors) &&
        all(factor -> all(isfinite, factor), decomposition.right_factors)
    diagnostics_are_finite =
        isfinite(coefficient_threshold) &&
        isfinite(normal_form_reconstruction_residual) &&
        isfinite(normal_form_reconstruction_tolerance) &&
        isfinite(filter_identity_residual) &&
        isfinite(filter_identity_tolerance)
    if !decomposition_is_finite ||
        !diagnostics_are_finite ||
        normal_form_reconstruction_residual > normal_form_reconstruction_tolerance ||
        filter_identity_residual > filter_identity_tolerance
        return _filter_normal_form_result(;
            filtered_operator=filtered_operator,
            centered_operator=centered_operator,
            coefficients=coefficients,
            left_operators=decomposition.left_factors,
            right_operators=decomposition.right_factors,
            local_filters=sinkhorn.local_filters,
            left_filter=sinkhorn.left_filter,
            right_filter=sinkhorn.right_filter,
            sinkhorn_result=sinkhorn,
            status=:numerical_failure,
            input_dimension=dimension,
            normal_form_trace=normal_form_trace,
            input_numerical_rank=input_numerical_rank,
            input_rank_threshold=input_rank_threshold,
            input_minimum_eigenvalue=minimum_input_eigenvalue,
            input_maximum_eigenvalue=maximum_input_eigenvalue,
            input_numerically_full_rank=input_numerical_rank == dimension,
            coefficient_numerical_rank=coefficient_numerical_rank,
            coefficient_threshold=coefficient_threshold,
            decomposition_hermiticity_residual=decomposition_hermiticity_residual,
            decomposition_hermiticity_tolerance=decomposition_hermiticity_tolerance,
            normal_form_reconstruction_residual=normal_form_reconstruction_residual,
            normal_form_reconstruction_tolerance=normal_form_reconstruction_tolerance,
            filter_identity_residual=filter_identity_residual,
            filter_identity_tolerance=filter_identity_tolerance,
            work_used=work_used,
            work_limit=budget.work_limit,
            estimated_postprocessing_work=budget.postprocessing_work,
            message="a checked normal-form reconstruction or local-filter identity " *
                    "residual exceeded its explicit verification tolerance",
        )
    end

    return _filter_normal_form_result(;
        filtered_operator=filtered_operator,
        centered_operator=centered_operator,
        coefficients=coefficients,
        left_operators=decomposition.left_factors,
        right_operators=decomposition.right_factors,
        local_filters=sinkhorn.local_filters,
        left_filter=sinkhorn.left_filter,
        right_filter=sinkhorn.right_filter,
        sinkhorn_result=sinkhorn,
        status=:converged,
        input_dimension=dimension,
        normal_form_trace=normal_form_trace,
        input_numerical_rank=input_numerical_rank,
        input_rank_threshold=input_rank_threshold,
        input_minimum_eigenvalue=minimum_input_eigenvalue,
        input_maximum_eigenvalue=maximum_input_eigenvalue,
        input_numerically_full_rank=input_numerical_rank == dimension,
        coefficient_numerical_rank=coefficient_numerical_rank,
        coefficient_threshold=coefficient_threshold,
        decomposition_hermiticity_residual=decomposition_hermiticity_residual,
        decomposition_hermiticity_tolerance=decomposition_hermiticity_tolerance,
        normal_form_reconstruction_residual=normal_form_reconstruction_residual,
        normal_form_reconstruction_tolerance=normal_form_reconstruction_tolerance,
        filter_identity_residual=filter_identity_residual,
        filter_identity_tolerance=filter_identity_tolerance,
        work_used=work_used,
        work_limit=budget.work_limit,
        estimated_postprocessing_work=budget.postprocessing_work,
        message="bounded Sinkhorn scaling, Hermitian operator-Schmidt decomposition, " *
                "and both reconstruction identities passed their checked tolerances",
    )
end

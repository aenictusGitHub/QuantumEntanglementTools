module QuantumEntanglementToolsEntanglementDetectionExt

using EntanglementDetection: EntanglementDetection
using LinearAlgebra: LinearAlgebra
using Random: Random
using Serialization: Serialization
import QuantumEntanglementTools as QET

const AUDITED_BACKEND_VERSION = v"0.2.2"
const WORKER_PATH = joinpath(@__DIR__, "entanglement_detection_worker.jl")
const PROCESS_OUTPUT_LIMIT = 8_192
const RESPONSE_BASE_LIMIT = 1_048_576
const RESPONSE_BYTES_PER_MATRIX_ENTRY = 64
const RESPONSE_ABSOLUTE_LIMIT = 268_435_456
const TERMINATION_GRACE_SECONDS = 0.1
const TERMINATION_FORCE_SECONDS = 2.0

backend_version() = Base.pkgversion(EntanglementDetection)

function _unknown_report(evidence, message::AbstractString)
    attempt = QET.EntanglementAttempt(
        :heuristic_search,
        :entanglement_detection,
        :unknown,
        false,
        nothing,
        evidence,
        message,
    )
    return QET.EntanglementReport(
        :unknown,
        false,
        nothing,
        :heuristic_search,
        :entanglement_detection,
        evidence,
        QET.EntanglementAttempt[attempt],
        message,
    )
end

function _failure_report(
    kind::Symbol,
    message::AbstractString;
    backend_message=nothing,
    backend_error_type=nothing,
    exitcode=nothing,
    termination=nothing,
    stdout="",
    stderr="",
)
    evidence = (
        schema_version=1,
        backend_version=backend_version(),
        backend_conclusion=:unavailable,
        certified_by_adapter=false,
        isolation=:child_process,
        execution=(
            ok=false,
            error_kind=kind,
            error_type=backend_error_type,
            backend_message=backend_message,
            exitcode=exitcode,
            termination=termination,
            captured_stdout=stdout,
            captured_stderr=stderr,
        ),
    )
    return _unknown_report(evidence, message)
end

function _output_excerpt(path::AbstractString)
    isfile(path) || return ""
    bytes = open(path, "r") do io
        return read(io, PROCESS_OUTPUT_LIMIT + 1)
    end
    length(bytes) <= PROCESS_OUTPUT_LIMIT && return String(bytes)
    return String(bytes[1:PROCESS_OUTPUT_LIMIT]) * "\n[output truncated by adapter]"
end

function _response_size_limit(request)
    entries = length(request.state)
    scalable_entries = fld(
        RESPONSE_ABSOLUTE_LIMIT - RESPONSE_BASE_LIMIT, RESPONSE_BYTES_PER_MATRIX_ENTRY
    )
    entries >= scalable_entries && return RESPONSE_ABSOLUTE_LIMIT
    return RESPONSE_BASE_LIMIT + RESPONSE_BYTES_PER_MATRIX_ENTRY * entries
end

function _decode_response(path::AbstractString, size_limit::Integer)
    response_size = filesize(path)
    response_size <= size_limit || return (
        ok=false,
        error_kind=:oversized_response,
        message="response size $response_size exceeds the adapter limit $size_limit",
        payload=nothing,
    )
    payload = try
        open(Serialization.deserialize, path)
    catch error
        return (
            ok=false,
            error_kind=:invalid_serialization,
            message=sprint(showerror, error),
            payload=nothing,
        )
    end
    return (ok=true, error_kind=nothing, message="", payload)
end

function _child_command(request_path::AbstractString, response_path::AbstractString)
    julia = Base.julia_cmd()
    active_project = Base.active_project()
    project_arguments =
        isnothing(active_project) ? String[] : ["--project=$(dirname(active_project))"]
    return `$julia --startup-file=no --history-file=no $project_arguments $WORKER_PATH $request_path $response_path`
end

function _process_has_exited(process)
    return try
        process_exited(process)
    catch
        false
    end
end

function _process_exited_within(process, timeout_seconds)
    _process_has_exited(process) && return true
    return try
        wait_status = timedwait(
            () -> _process_has_exited(process),
            timeout_seconds;
            pollint=max(0.001, min(0.01, timeout_seconds / 10)),
        )
        wait_status === :ok
    catch
        _process_has_exited(process)
    end
end

function _terminate_child(process)
    graceful_error = nothing
    graceful_signal_sent = try
        kill(process, Base.SIGTERM)
        true
    catch error
        graceful_error = sprint(showerror, error)
        false
    end
    if _process_exited_within(process, TERMINATION_GRACE_SECONDS)
        return (;
            graceful_signal_sent,
            graceful_error,
            force_signal_sent=false,
            force_signal=nothing,
            force_error=nothing,
            reaped=true,
        )
    end

    # Libuv maps SIGTERM to TerminateProcess on Windows. On Unix, escalate a
    # TERM-resistant child to the uncatchable SIGKILL.
    force_signal = Sys.iswindows() ? Base.SIGTERM : Base.SIGKILL
    force_error = nothing
    force_signal_sent = try
        kill(process, force_signal)
        true
    catch error
        force_error = sprint(showerror, error)
        false
    end
    reaped = _process_exited_within(process, TERMINATION_FORCE_SECONDS)
    return (;
        graceful_signal_sent,
        graceful_error,
        force_signal_sent,
        force_signal,
        force_error,
        reaped,
    )
end

function _wait_for_process(
    process,
    timeout_seconds;
    child_state::Base.RefValue{Symbol}=Ref(:running),
    wait_function=timedwait,
)
    child_state[] = :running
    lifecycle_handled = false
    try
        wait_status = wait_function(
            () -> process_exited(process),
            timeout_seconds;
            pollint=max(0.001, min(0.05, timeout_seconds / 10)),
        )
        wait_status in (:ok, :timed_out) ||
            throw(ErrorException("process wait returned unexpected status $wait_status"))
        if wait_status === :timed_out
            termination = _terminate_child(process)
            child_state[] = termination.reaped ? :reaped : :unreaped
            lifecycle_handled = true
            return (
                status=termination.reaped ? :timed_out : :termination_failed,
                process,
                termination,
                error_type=nothing,
                message=if termination.reaped
                    "the isolated backend process exceeded its wall-clock limit"
                else
                    "the isolated backend process exceeded its wall-clock limit and could not be reaped after forced termination"
                end,
            )
        end

        _process_has_exited(process) ||
            throw(ErrorException("process wait completed before the child exited"))
        child_state[] = :reaped
        lifecycle_handled = true
        return (
            status=success(process) ? :ok : :process_failed,
            process,
            termination=nothing,
            error_type=nothing,
            message=if success(process)
                ""
            else
                "the isolated backend process exited unsuccessfully"
            end,
        )
    finally
        if !lifecycle_handled
            termination = _terminate_child(process)
            child_state[] = termination.reaped ? :reaped : :unreaped
        end
    end
end

function _wait_for_child(
    command::Cmd,
    stdout_path,
    stderr_path,
    timeout_seconds;
    child_state::Base.RefValue{Symbol}=Ref(:not_started),
    wait_function=timedwait,
    close_function=close,
)
    child_stdout = nothing
    child_stderr = nothing
    process = nothing
    setup_error = nothing
    try
        child_stdout = open(stdout_path, "w")
        child_stderr = open(stderr_path, "w")
        process = run(
            pipeline(command; stdout=child_stdout, stderr=child_stderr); wait=false
        )
        child_state[] = :running
    catch error
        setup_error = error
    finally
        for stream in (child_stderr, child_stdout)
            isnothing(stream) && continue
            try
                close_function(stream)
            catch error
                isnothing(setup_error) && (setup_error = error)
            end
        end
    end

    if !isnothing(setup_error)
        termination = if isnothing(process)
            child_state[] = :not_started
            nothing
        else
            result = _terminate_child(process)
            child_state[] = result.reaped ? :reaped : :unreaped
            result
        end
        return (
            status=:launch_failed,
            process,
            termination,
            error_type=string(typeof(setup_error)),
            message=sprint(showerror, setup_error),
        )
    end

    return _wait_for_process(process, timeout_seconds; child_state, wait_function)
end

function _validated_request(
    rho::AbstractMatrix{<:Number}, dims, method::QET.EntanglementDetectionSearch
)
    backend_version() == AUDITED_BACKEND_VERSION || throw(
        ArgumentError(
            "the adapter was audited only for EntanglementDetection.jl " *
            "$AUDITED_BACKEND_VERSION; loaded version is $(backend_version())",
        ),
    )
    layout = QET.SubsystemLayout(dims)
    length(layout) >= 2 ||
        throw(ArgumentError("EntanglementDetectionSearch requires at least two subsystems"))
    size(rho) == (layout.total_dimension, layout.total_dimension) || throw(
        DimensionMismatch(
            "prod(dims)=$(layout.total_dimension) does not match state size $(size(rho))",
        ),
    )
    analysis = QET._tierd_density_analysis(
        rho;
        atol=method.atol,
        rtol=method.rtol,
        allow_densify=method.allow_densify,
        boundary_policy=:reject,
        operation="EntanglementDetectionSearch",
    )
    return (
        schema_version=1,
        backend_version=AUDITED_BACKEND_VERSION,
        state=analysis.matrix,
        dims=Tuple(layout),
        max_iteration=method.max_iteration,
        epsilon=method.epsilon,
        callback_iter=method.callback_iter,
    )
end

function _invalid_response(message::AbstractString, stdout::String, stderr::String)
    return _failure_report(
        :invalid_response,
        "the isolated backend returned an invalid response";
        backend_message=message,
        stdout,
        stderr,
    )
end

_has_payload_fields(payload, fields) = all(field -> hasproperty(payload, field), fields)

function _translate_response(
    payload,
    stdout::String,
    stderr::String,
    expected_size::Tuple{Int,Int},
    expected_representation::Type,
)
    payload isa NamedTuple ||
        return _invalid_response("response is not a NamedTuple", stdout, stderr)
    common_fields = (:schema_version, :backend_version, :ok)
    _has_payload_fields(payload, common_fields) ||
        return _invalid_response("response omitted common metadata", stdout, stderr)
    payload.schema_version isa Integer &&
    !(payload.schema_version isa Bool) &&
    payload.schema_version == 1 ||
        return _invalid_response("response schema version is not 1", stdout, stderr)
    payload.backend_version isa VersionNumber &&
    payload.backend_version == AUDITED_BACKEND_VERSION ||
        return _invalid_response("response backend version is invalid", stdout, stderr)
    payload.ok isa Bool ||
        return _invalid_response("response status is not Boolean", stdout, stderr)

    if !payload.ok
        error_fields = (:error_type, :message)
        _has_payload_fields(payload, error_fields) ||
            return _invalid_response("error response omitted diagnostics", stdout, stderr)
        payload.error_type isa AbstractString ||
            return _invalid_response("error_type is not text", stdout, stderr)
        payload.message isa AbstractString ||
            return _invalid_response("error message is not text", stdout, stderr)
        return _failure_report(
            :backend_exception,
            "EntanglementDetection.jl raised an exception; no mathematical conclusion was made";
            backend_message=String(payload.message),
            backend_error_type=String(payload.error_type),
            stdout,
            stderr,
        )
    end

    success_fields = (
        :backend_conclusion,
        :backend_input_representation,
        :witness_operator,
        :witness_expectation,
        :witness_hermiticity_residual,
        :active_set_size,
        :rng_state_changed,
        :stdout_binding_changed,
        :logger_changed,
        :blas_threads_before,
        :blas_threads_after,
    )
    _has_payload_fields(payload, success_fields) ||
        return _invalid_response("success response omitted evidence", stdout, stderr)
    candidate = payload.backend_conclusion
    candidate isa Symbol && candidate in (:entangled, :separable, :inconclusive) ||
        return _invalid_response("candidate conclusion is invalid", stdout, stderr)
    payload.backend_input_representation === expected_representation ||
        return _invalid_response("backend input representation is invalid", stdout, stderr)

    witness = payload.witness_operator
    witness isa AbstractMatrix{<:Number} && size(witness) == expected_size ||
        return _invalid_response(
            "witness operator shape or element type is invalid", stdout, stderr
        )
    all(isfinite, witness) || return _invalid_response(
        "witness operator contains non-finite entries", stdout, stderr
    )
    expectation = payload.witness_expectation
    expectation isa Real && !(expectation isa Bool) && isfinite(expectation) ||
        return _invalid_response("witness expectation is invalid", stdout, stderr)
    hermiticity_residual = payload.witness_hermiticity_residual
    hermiticity_residual isa Real &&
    !(hermiticity_residual isa Bool) &&
    isfinite(hermiticity_residual) &&
    hermiticity_residual >= zero(hermiticity_residual) ||
        return _invalid_response("witness Hermiticity residual is invalid", stdout, stderr)
    active_set_size = payload.active_set_size
    active_set_size isa Integer && !(active_set_size isa Bool) && active_set_size >= 0 ||
        return _invalid_response("active-set size is invalid", stdout, stderr)

    state_flags = (
        payload.rng_state_changed, payload.stdout_binding_changed, payload.logger_changed
    )
    all(flag -> flag isa Bool, state_flags) ||
        return _invalid_response("child-state flags are invalid", stdout, stderr)
    blas_threads = (payload.blas_threads_before, payload.blas_threads_after)
    all(
        threads -> threads isa Integer && !(threads isa Bool) && threads > 0, blas_threads
    ) || return _invalid_response("child BLAS thread counts are invalid", stdout, stderr)

    evidence = (
        schema_version=1,
        backend_version=payload.backend_version,
        backend_conclusion=candidate,
        certified_by_adapter=false,
        witness=(operator=witness, expectation, hermiticity_residual),
        active_set_size,
        backend_input_representation=payload.backend_input_representation,
        isolation=:child_process,
        execution=(
            ok=true,
            exitcode=0,
            captured_stdout=stdout,
            captured_stderr=stderr,
            child_rng_state_changed=payload.rng_state_changed,
            child_stdout_binding_changed=payload.stdout_binding_changed,
            child_logger_changed=payload.logger_changed,
            child_blas_threads_before=payload.blas_threads_before,
            child_blas_threads_after=payload.blas_threads_after,
        ),
    )
    message = if candidate === :entangled
        "EntanglementDetection.jl produced an entangled candidate in an isolated process; the adapter did not independently validate a witness certificate"
    elseif candidate === :separable
        "EntanglementDetection.jl produced a separable candidate in an isolated process; the adapter did not independently validate a decomposition certificate"
    else
        "EntanglementDetection.jl was inconclusive in an isolated process"
    end
    return _unknown_report(evidence, message)
end

function _detect_entanglement_isolated(request, method::QET.EntanglementDetectionSearch)
    directory = mktempdir(; prefix="qet-entanglement-detection-")
    cleanup_synchronously = true
    child_state = Ref(:not_started)
    try
        request_path = joinpath(directory, "request.bin")
        response_path = joinpath(directory, "response.bin")
        stdout_path = joinpath(directory, "stdout.log")
        stderr_path = joinpath(directory, "stderr.log")
        response_limit = _response_size_limit(request)
        open(request_path, "w") do io
            return Serialization.serialize(io, request)
        end

        process_result = _wait_for_child(
            _child_command(request_path, response_path),
            stdout_path,
            stderr_path,
            method.timeout_seconds,
            ;
            child_state,
        )
        captured_stdout = _output_excerpt(stdout_path)
        captured_stderr = _output_excerpt(stderr_path)
        if process_result.status === :launch_failed
            return _failure_report(
                :launch_failed,
                "the isolated EntanglementDetection.jl process or its output capture could not be initialized";
                backend_message=process_result.message,
                backend_error_type=process_result.error_type,
                termination=process_result.termination,
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        elseif process_result.status === :timed_out
            return _failure_report(
                :timeout,
                "the isolated EntanglementDetection.jl process timed out; no mathematical conclusion was made";
                termination=process_result.termination,
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        elseif process_result.status === :termination_failed
            cleanup_synchronously = false
            return _failure_report(
                :termination_failed,
                "the isolated EntanglementDetection.jl process timed out and resisted bounded termination; no mathematical conclusion was made";
                termination=process_result.termination,
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        elseif process_result.status === :process_failed
            return _failure_report(
                :process_failed,
                "the isolated EntanglementDetection.jl process failed; no mathematical conclusion was made";
                exitcode=process_result.process.exitcode,
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        elseif !isfile(response_path)
            return _failure_report(
                :missing_response,
                "the isolated EntanglementDetection.jl process returned no response";
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        end

        decoded = _decode_response(response_path, response_limit)
        if !decoded.ok
            return _failure_report(
                :invalid_response,
                "the isolated EntanglementDetection.jl response could not be decoded";
                backend_message=decoded.message,
                backend_error_type=String(decoded.error_kind),
                stdout=captured_stdout,
                stderr=captured_stderr,
            )
        end
        expected_representation = typeof(complex(zero(eltype(request.state))))
        return _translate_response(
            decoded.payload,
            captured_stdout,
            captured_stderr,
            size(request.state),
            expected_representation,
        )
    finally
        child_state[] === :unreaped && (cleanup_synchronously = false)
        if cleanup_synchronously && ispath(directory)
            try
                rm(directory; recursive=true)
            catch
                # `mktempdir` registered a process-exit cleanup as a fallback.
            end
        end
    end
end

function QET.detect_entanglement(
    rho::AbstractMatrix{<:Number}, dims, method::QET.EntanglementDetectionSearch
)
    request = _validated_request(rho, dims, method)
    caller_rng = copy(Random.default_rng())
    try
        return _detect_entanglement_isolated(request, method)
    finally
        copy!(Random.default_rng(), caller_rng)
    end
end

end

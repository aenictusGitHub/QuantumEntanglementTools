#!/usr/bin/env julia

using LinearAlgebra
using Logging
using QuantumEntanglementTools
using Random
using Serialization
using SparseArrays
using Test
using EntanglementDetection

const QETED = QuantumEntanglementTools
const LOAD_ORDER_PROBE = joinpath(@__DIR__, "load_order_probe.jl")

function extension_module()
    extension = Base.get_extension(QETED, :QuantumEntanglementToolsEntanglementDetectionExt)
    isnothing(extension) && error("EntanglementDetection extension is not loaded")
    return extension
end

function run_load_order_probe(order)
    project_file = Base.active_project()
    isnothing(project_file) && error("no active test project")
    command = `$(Base.julia_cmd()) --startup-file=no --history-file=no --project=$(dirname(project_file)) $LOAD_ORDER_PROBE $order`
    return success(pipeline(command; stdout=devnull, stderr=stderr))
end

@testset "EntanglementDetection optional extension" begin
    @testset "load order and capabilities" begin
        @test run_load_order_probe("core_first")
        @test run_load_order_probe("backend_first")

        extension = extension_module()
        backends = QETED.available_entanglement_backends()
        @test length(backends) == 2
        @test first(backends) isa QETED.NativeEntanglementBackend
        @test last(backends) isa QETED.EntanglementDetectionBackend
        capabilities = QETED.backend_capabilities(last(backends))
        @test capabilities.version == v"0.2.2"
        @test capabilities.methods == (:heuristic_search,)
        @test capabilities.conclusions == (:unknown,)
        @test capabilities.backend_candidates == (:entangled, :separable, :inconclusive)
        @test !capabilities.side_effect_free
        @test capabilities.caller_state_isolated
        @test capabilities.isolation === :child_process
        @test capabilities.minimum_resolvable_julia == v"1.11.0"
        @test !capabilities.certifies_conclusions &&
            !capabilities.source_integrity_enforced &&
            capabilities.transport_trust === :same_version_local_worker &&
            !capabilities.resource_sandboxed
        status = QETED.backend_status()
        optional = status.entanglement.entanglement_detection
        @test optional.loaded
        @test optional.ready
        @test optional.version == v"0.2.2"
        @test optional.dependency === :EntanglementDetection
        @test occursin("uncertified", optional.message)
        @test isempty(
            Test.detect_ambiguities(
                QETED, EntanglementDetection, extension; recursive=false
            ),
        )
    end

    @testset "configuration and input validation" begin
        method = QETED.EntanglementDetectionSearch(
            timeout_seconds=30,
            max_iteration=2,
            epsilon=1.0f-5,
            callback_iter=2,
            atol=0,
            rtol=1.0f-6,
            allow_densify=true,
        )
        @test method.timeout_seconds == 30
        @test method.max_iteration == 2
        @test method.epsilon === 1.0f-5
        @test method.callback_iter == 2
        @test method.atol === 0.0f0
        @test method.rtol === 1.0f-6
        @test method.allow_densify

        @test_throws ArgumentError QETED.EntanglementDetectionSearch(timeout_seconds=0)
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(timeout_seconds=Inf)
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(
            timeout_seconds=big"1e-1000"
        )
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(max_iteration=0)
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(callback_iter=false)
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(epsilon=-1)
        @test_throws ArgumentError QETED.EntanglementDetectionSearch(atol=NaN)

        mixed = Matrix{Float64}(I, 4, 4) / 4
        rng_before_validation = copy(Random.default_rng())
        @test_throws ArgumentError QETED.detect_entanglement(
            sparse(mixed), (2, 2), QETED.EntanglementDetectionSearch()
        )
        @test copy(Random.default_rng()) == rng_before_validation
        @test_throws DimensionMismatch QETED.detect_entanglement(
            mixed, (2, 3), QETED.EntanglementDetectionSearch()
        )
        @test_throws ArgumentError QETED.detect_entanglement(
            mixed, (4,), QETED.EntanglementDetectionSearch()
        )
        @test_throws ArgumentError QETED.detect_entanglement(
            2mixed, (2, 2), QETED.EntanglementDetectionSearch()
        )
    end

    @testset "response validation and bounded IPC" begin
        extension = extension_module()
        expected_size = (4, 4)
        expected_representation = ComplexF64
        invalid_report(payload) = extension._translate_response(
            payload, "", "", expected_size, expected_representation
        )
        common = (schema_version=1, backend_version=v"0.2.2")

        malformed_payloads = (
            merge(common, (ok=missing,)),
            merge(common, (ok=1,)),
            (schema_version=missing, backend_version=v"0.2.2", ok=false),
            (schema_version=true, backend_version=v"0.2.2", ok=false),
            (schema_version=2, backend_version=v"0.2.2", ok=false),
            (schema_version=1, backend_version=v"0.2.1", ok=false),
            merge(common, (ok=false, error_type=1, message="failure")),
            merge(common, (ok=false, error_type="ArgumentError", message=[])),
            merge(common, (ok=true,)),
        )
        for payload in malformed_payloads
            report = invalid_report(payload)
            @test report.status === :unknown
            @test !report.certified
            @test report.evidence.execution.error_kind === :invalid_response
        end

        valid_success = merge(
            common,
            (
                ok=true,
                backend_conclusion=:inconclusive,
                backend_input_representation=ComplexF64,
                witness_operator=Matrix{ComplexF64}(I, 4, 4),
                witness_expectation=0.0,
                witness_hermiticity_residual=0.0,
                active_set_size=1,
                rng_state_changed=true,
                stdout_binding_changed=false,
                logger_changed=false,
                blas_threads_before=1,
                blas_threads_after=1,
            ),
        )
        for payload in (
            merge(valid_success, (backend_conclusion=:unsupported,)),
            merge(valid_success, (witness_operator="not a matrix",)),
            merge(valid_success, (witness_expectation=true,)),
            merge(valid_success, (witness_hermiticity_residual=true,)),
            merge(valid_success, (rng_state_changed=1,)),
            merge(valid_success, (active_set_size=-1,)),
        )
            report = invalid_report(payload)
            @test report.evidence.execution.error_kind === :invalid_response
        end

        backend_error = invalid_report(
            merge(
                common,
                (ok=false, error_type="ArgumentError", message="injected backend failure"),
            ),
        )
        @test backend_error.status === :unknown
        @test backend_error.evidence.execution.error_kind === :backend_exception

        request = (state=zeros(ComplexF64, 4, 4),)
        response_limit = extension._response_size_limit(request)
        @test response_limit ==
            extension.RESPONSE_BASE_LIMIT + extension.RESPONSE_BYTES_PER_MATRIX_ENTRY * 16
        mktempdir() do directory
            oversized_path = joinpath(directory, "oversized.bin")
            open(oversized_path, "w") do io
                write(io, zeros(UInt8, 33))
            end
            oversized = extension._decode_response(oversized_path, 32)
            @test !oversized.ok
            @test oversized.error_kind === :oversized_response

            malformed_path = joinpath(directory, "malformed.bin")
            serialized_buffer = IOBuffer()
            serialize(serialized_buffer, fill(UInt8(0x42), 1_024))
            serialized_bytes = take!(serialized_buffer)
            open(malformed_path, "w") do io
                write(io, serialized_bytes[1:(end - 100)])
            end
            malformed = extension._decode_response(malformed_path, length(serialized_bytes))
            @test !malformed.ok
            @test malformed.error_kind === :invalid_serialization

            valid_path = joinpath(directory, "valid.bin")
            open(valid_path, "w") do io
                serialize(io, valid_success)
            end
            decoded = extension._decode_response(valid_path, response_limit)
            @test decoded.ok
            @test decoded.payload == valid_success

            output_path = joinpath(directory, "output.log")
            open(output_path, "w") do io
                write(io, fill(UInt8('x'), extension.PROCESS_OUTPUT_LIMIT + 100))
            end
            excerpt = extension._output_excerpt(output_path)
            @test occursin("[output truncated by adapter]", excerpt)
            @test ncodeunits(excerpt) < extension.PROCESS_OUTPUT_LIMIT + 100
        end
    end

    @testset "bounded process lifecycle" begin
        extension = extension_module()

        for injected_error in (InterruptException(), ErrorException("injected wait error"))
            process = run(
                `$(Base.julia_cmd()) --startup-file=no --history-file=no -e $("sleep(30)")`;
                wait=false,
            )
            child_state = Ref(:running)
            wait_function = (callback, timeout; pollint) -> throw(injected_error)
            @test_throws typeof(injected_error) extension._wait_for_process(
                process, 30.0; child_state, wait_function
            )
            @test child_state[] === :reaped
            @test process_exited(process)
        end

        mktempdir() do directory
            child_state = Ref(:not_started)
            close_count = Ref(0)
            close_after_throwing = function (stream)
                close(stream)
                close_count[] += 1
                close_count[] == 1 &&
                    throw(ErrorException("injected post-launch stream-close failure"))
                return nothing
            end
            result = extension._wait_for_child(
                `$(Base.julia_cmd()) --startup-file=no --history-file=no -e $("sleep(30)")`,
                joinpath(directory, "stdout.log"),
                joinpath(directory, "stderr.log"),
                30.0;
                child_state,
                close_function=close_after_throwing,
            )
            @test result.status === :launch_failed &&
                !isnothing(result.process) &&
                result.termination.reaped &&
                child_state[] === :reaped &&
                process_exited(result.process)
        end

        if Sys.isunix()
            directory_path = Ref("")
            termination, elapsed, process = mktempdir() do directory
                directory_path[] = directory
                ready_path = joinpath(directory, "ready")
                script = "trap '' TERM; printf ready > \"\$QET_READY\"; exec sleep 30"
                command = setenv(Cmd(["/bin/sh", "-c", script]), "QET_READY" => ready_path)
                child = run(command; wait=false)
                ready_status = timedwait(() -> isfile(ready_path), 2.0; pollint=0.01)
                @test ready_status === :ok
                started = time()
                result = extension._terminate_child(child)
                return result, time() - started, child
            end
            @test termination.graceful_signal_sent
            @test termination.force_signal_sent
            @test termination.force_signal == Base.SIGKILL
            @test termination.reaped
            @test process_exited(process)
            @test elapsed <
                extension.TERMINATION_GRACE_SECONDS +
                  extension.TERMINATION_FORCE_SECONDS +
                  1
            @test !ispath(directory_path[])
        end
    end

    @testset "isolated search and caller-state preservation" begin
        mixed = Matrix{Float64}(I, 4, 4) / 4
        method = QETED.EntanglementDetectionSearch(
            timeout_seconds=120, max_iteration=2, epsilon=1e-5, callback_iter=2
        )

        rng_snapshot = copy(Random.default_rng())
        expected_rng = copy(rng_snapshot)
        expected_values = rand(expected_rng, UInt64, 8)
        stdout_before = stdout
        stderr_before = stderr
        logger_before = Logging.global_logger()
        blas_threads_before = BLAS.get_num_threads()

        report = QETED.detect_entanglement(mixed, (2, 2), method)

        @test report.status === :unknown
        @test !report.certified
        @test report.certificate_kind === nothing
        @test report.method === :heuristic_search
        @test report.backend === :entanglement_detection
        @test length(report.attempts) == 1
        @test !only(report.attempts).certified
        @test report.evidence.certified_by_adapter === false
        @test report.evidence.isolation === :child_process
        @test report.evidence.execution.ok
        @test report.evidence.backend_conclusion in (:entangled, :separable, :inconclusive)
        @test report.evidence.backend_input_representation === ComplexF64
        @test report.evidence.witness.operator isa Matrix
        @test isfinite(report.evidence.witness.expectation)
        @test report.evidence.execution.child_rng_state_changed === true
        @test copy(Random.default_rng()) == rng_snapshot
        @test rand(Random.default_rng(), UInt64, 8) == expected_values
        @test stdout === stdout_before
        @test stderr === stderr_before
        @test Logging.global_logger() === logger_before
        @test BLAS.get_num_threads() == blas_threads_before
    end

    @testset "timeout is an unknown execution result" begin
        mixed = Matrix{Float64}(I, 4, 4) / 4
        rng_before = copy(Random.default_rng())
        stdout_before = stdout
        logger_before = Logging.global_logger()
        blas_threads_before = BLAS.get_num_threads()
        started = time()
        report = QETED.detect_entanglement(
            mixed,
            (2, 2),
            QETED.EntanglementDetectionSearch(
                timeout_seconds=1e-6, max_iteration=1, callback_iter=1
            ),
        )
        elapsed = time() - started
        @test report.status === :unknown
        @test !report.certified
        @test report.evidence.backend_conclusion === :unavailable
        @test !report.evidence.execution.ok
        @test report.evidence.execution.error_kind === :timeout
        @test report.evidence.execution.termination.reaped
        extension = extension_module()
        @test elapsed <
            1e-6 +
              extension.TERMINATION_GRACE_SECONDS +
              extension.TERMINATION_FORCE_SECONDS +
              1
        @test copy(Random.default_rng()) == rng_before
        @test stdout === stdout_before
        @test Logging.global_logger() === logger_before
        @test BLAS.get_num_threads() == blas_threads_before
    end
end

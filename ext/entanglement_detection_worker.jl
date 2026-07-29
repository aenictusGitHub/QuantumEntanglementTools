#!/usr/bin/env julia

using EntanglementDetection: EntanglementDetection
using LinearAlgebra: LinearAlgebra
using Logging: Logging
using Random: Random
using Serialization: Serialization

const AUDITED_BACKEND_VERSION = v"0.2.2"

function backend_response(request)
    request.schema_version == 1 ||
        throw(ArgumentError("unsupported adapter request schema"))
    Base.pkgversion(EntanglementDetection) == request.backend_version || throw(
        ArgumentError(
            "backend version mismatch: expected $(request.backend_version), " *
            "loaded $(Base.pkgversion(EntanglementDetection))",
        ),
    )
    request.backend_version == AUDITED_BACKEND_VERSION ||
        throw(ArgumentError("worker is restricted to backend version 0.2.2"))

    rng_before = copy(Random.default_rng())
    stdout_before = stdout
    logger_before = Logging.global_logger()
    blas_threads_before = LinearAlgebra.BLAS.get_num_threads()

    # Version 0.2.2 produces a complex separable estimate even for real input,
    # then requires that estimate and the input to share an element type.
    # Convert representation without changing the real component precision.
    backend_state = complex.(request.state)
    result = EntanglementDetection.entanglement_detection(
        backend_state,
        request.dims;
        measure="2-norm",
        max_iteration=request.max_iteration,
        epsilon=request.epsilon,
        callback_iter=request.callback_iter,
        verbose=0,
        logfile=nothing,
        shortcut=false,
    )

    witness_operator = Matrix(result.witness.W)
    witness_expectation = real(LinearAlgebra.dot(witness_operator, backend_state))
    hermiticity_residual = maximum(
        abs, witness_operator - adjoint(witness_operator); init=0.0
    )
    conclusion = if result.ent === true
        :entangled
    elseif result.ent === false
        :separable
    else
        :inconclusive
    end
    return (
        schema_version=1,
        ok=true,
        backend_version=Base.pkgversion(EntanglementDetection),
        backend_conclusion=conclusion,
        backend_input_representation=eltype(backend_state),
        witness_operator=witness_operator,
        witness_expectation=witness_expectation,
        witness_hermiticity_residual=hermiticity_residual,
        active_set_size=length(result.decompose),
        rng_state_changed=copy(Random.default_rng()) != rng_before,
        stdout_binding_changed=stdout !== stdout_before,
        logger_changed=Logging.global_logger() !== logger_before,
        blas_threads_before,
        blas_threads_after=LinearAlgebra.BLAS.get_num_threads(),
    )
end

function exception_response(error)
    return (
        schema_version=1,
        ok=false,
        backend_version=Base.pkgversion(EntanglementDetection),
        error_type=string(typeof(error)),
        message=sprint(showerror, error),
    )
end

function main(arguments)
    length(arguments) == 2 ||
        throw(ArgumentError("expected request and response file paths"))
    request_path, response_path = arguments
    response = try
        request = open(Serialization.deserialize, request_path)
        backend_response(request)
    catch error
        exception_response(error)
    end
    open(response_path, "w") do io
        return Serialization.serialize(io, response)
    end
    return nothing
end

main(ARGS)

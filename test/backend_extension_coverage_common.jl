using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using Random

function run_backend_extension_coverage(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
    Backends.disable_scalar_backend!(B)
    ArrayBackend = Backends.gpu_backend_array_type(B)

    @test Backends.gpu_backend_loaded(B)
    @test ArrayBackend !== nothing
    @test Backends.gpu_backend_name(B) !== nothing
    @test Backends.gpu_backend_name(ArrayBackend) !== nothing
    @test Backends.array_backend_selector(
        ArrayBackend) isa Backends.AbstractArrayBackend
    @test B in Backends.available_gpu_backends()

    uniform = Backends.backend_rand(B, Float32, 2, 2)
    normal = Backends.backend_randn(B, Float32, 2, 2)
    zeros_array = Backends.backend_zeros(B, Float32, 2, 2)
    filled = Backends.backend_fill(B, 3.0f0, 2, 2)
    @test uniform isa ArrayBackend
    @test normal isa ArrayBackend
    @test zeros_array isa ArrayBackend
    @test filled isa ArrayBackend
    @test all(isfinite, Array(uniform))
    @test all(isfinite, Array(normal))
    @test iszero(sum(Array(zeros_array)))
    @test Array(filled) == fill(3.0f0, 2, 2)
    @test !isnothing(Backends.compute_device_identifier(zeros_array))

    context = Backends._prepare_device_execution_context(zeros_array)
    @test Backends._prepared_device_execution_compute_device(context) ==
        Backends.compute_device(zeros_array)
    result = Backends._with_prepared_device_execution_context(
        context,
    ) do
        fill!(zeros_array, 2.0f0)
        zeros_array
    end
    @test result === zeros_array
    Backends._synchronize_prepared_device_execution_context!(context)
    @test Array(zeros_array) == fill(2.0f0, 2, 2)

    basis_host = reshape(Float32.(1:12), 2, 3, 2)
    pupil_support_host = Bool[true false true; true true false]
    coefficients_host = Float32[0.25, -0.5]
    expansion = AdaptiveOpticsSim.Calibration.ModalOPDExpansionPlan(
        ArrayBackend(basis_host),
        ArrayBackend(pupil_support_host),
    )
    coefficients = ArrayBackend(coefficients_host)
    opd = ArrayBackend(zeros(Float32, 2, 3))
    @test AdaptiveOpticsSim.Calibration.combine_basis!(
        opd,
        expansion,
        coefficients,
    ) === opd
    Backends.synchronize_backend!(Backends.execution_style(opd))
    @test Array(opd) == Float32[-3.25 0.0 -4.25; -3.5 -4.0 0.0]
    return nothing
end

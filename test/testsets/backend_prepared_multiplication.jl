function prepared_mul_fixture(::Type{T}, shape) where {T}
    return reshape(T[convert(T, sin(0.41i)) +
        (T <: Complex ? convert(T, 0.3im * cos(0.27i)) : zero(T))
        for i in 1:prod(shape)], shape)
end

prepared_mul_allocations(C, A, B, coeffs) =
    @allocated Backends._enqueue_mul!(C, A, B, coeffs)

@testset "Prepared backend multiplication" begin
    for T in (Float16, Float32, Float64, ComplexF32, ComplexF64)
        A = prepared_mul_fixture(T, (13, 11))
        right = prepared_mul_fixture(T, (11, 9))
        transposed = transpose(prepared_mul_fixture(T, (9, 11)))
        vector = vec(prepared_mul_fixture(T, (11,)))
        for B in (right, transposed, vector)
            expected = A * B
            C = similar(expected)
            coeffs = Backends._prepare_mul_coefficients(C)
            @test coeffs === (; alpha=one(T), beta=zero(T))
            fill!(C, T(NaN))
            @test Backends._enqueue_mul!(C, A, B, coeffs) === C
            @test C ≈ expected rtol=32eps(real(T)) atol=32eps(real(T))
            Backends._enqueue_mul!(C, A, B, coeffs)
            prepared_mul_allocations(C, A, B, coeffs)
            if !coverage_instrumented()
                @test prepared_mul_allocations(C, A, B, coeffs) == 0
            end
        end
    end
end

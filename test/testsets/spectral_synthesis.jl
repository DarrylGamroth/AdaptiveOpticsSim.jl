import FixedSizeArrays

struct SpectralSynthesisForeignMatrix{T} <: AbstractMatrix{T}
    values::Matrix{T}
end
Base.size(matrix::SpectralSynthesisForeignMatrix) = size(matrix.values)
Base.getindex(matrix::SpectralSynthesisForeignMatrix, i::Int, j::Int) = matrix.values[i, j]
AdaptiveOpticsSim.Backends.backend(::SpectralSynthesisForeignMatrix) =
    AdaptiveOpticsSim.Backends.CUDABackend()

struct SpectralSynthesisShiftedMatrix{T} <: AbstractMatrix{T}
    values::Matrix{T}
end
Base.size(matrix::SpectralSynthesisShiftedMatrix) = size(matrix.values)
Base.axes(matrix::SpectralSynthesisShiftedMatrix) =
    (0:(size(matrix, 1) - 1), 0:(size(matrix, 2) - 1))
Base.getindex(matrix::SpectralSynthesisShiftedMatrix, i::Int, j::Int) =
    matrix.values[i + 1, j + 1]

function spectral_synthesis_allocations(owner)
    Optics.synthesize_spectrum!(owner)
    return @allocated Optics.synthesize_spectrum!(owner)
end

# Independent dense scalar calculation, rather than the execution loop's
# diagonal/pair decomposition, also exercises the complex coherence convention.
function spectral_synthesis_oracle(envelopes, basis, coefficients, visibility,
        channels, response)
    n_anchor = size(basis, 2)
    n_path = size(coefficients, 2)
    output = zeros(real(eltype(first(envelopes))), size(first(envelopes)))
    for index in eachindex(output), c in eachindex(channels.samples)
        fields = [coefficients[c, k] * sum(basis[c, a] *
            envelopes[(k - 1) * n_anchor + a][index] for a in 1:n_anchor)
            for k in 1:n_path]
        intensity = sum(fields[k] * visibility[c][k, l] * conj(fields[l])
            for k in 1:n_path, l in 1:n_path)
        output[index] += real(intensity) * channels[c].weight * response[c]
    end
    return output
end

@testset "Prepared spectral synthesis" begin
    for T in (Float32, Float64)
        channels = SpectralBundle(T[600e-9, 800e-9], T[1, 3]; T=T)
        envelopes = [fill(Complex{T}(2, 1), 2, 3),
            fill(Complex{T}(4, -1), 2, 3)]
        basis = T[1 0; 0 1]
        coefficients = reshape(Complex{T}[1, 2im], 2, 1)
        visibility = [ones(Complex{T}, 1, 1) for _ in 1:2]
        response = T[0.5, 0.25]
        expected = spectral_synthesis_oracle(envelopes, basis, coefficients,
            visibility, channels, response)
        owner = Optics.prepare_spectral_synthesis(envelopes, basis,
            coefficients, visibility, channels; response=response)
        @test @inferred(Optics.synthesize_spectrum!(owner)) === owner.products.output
        @test owner.products.output ≈ expected
        @test owner.plan.basis isa FixedSizeArrays.FixedSizeMatrix
        @test owner.plan.channels isa FixedSizeVector
        @test isconcretetype(eltype(owner.envelopes))
        @test isconcretetype(eltype(owner.workspace.channel_fields))
        @test isconcretetype(eltype(owner.plan.visibility))
        @test eltype(owner.products.output) === T
        # Calibration snapshots stay fixed; exact envelopes remain live bindings.
        fill!(basis, zero(T))
        fill!(coefficients, zero(Complex{T}))
        fill!(response, zero(T))
        fill!(visibility[1], zero(Complex{T}))
        channels.samples[1] = SpectralSample{T}(T(600e-9), zero(T))
        @test Optics.synthesize_spectrum!(owner) ≈ expected
        envelopes[1] .*= T(2)
        @test Optics.synthesize_spectrum!(owner) != expected
        bound = owner.envelopes[1]
        envelopes[1] = zeros(Complex{T}, 2, 3)
        @test owner.envelopes[1] === bound
        if !coverage_instrumented()
            @test spectral_synthesis_allocations(owner) == 0
        end

        channel = SpectralBundle(T[700e-9], T[1]; T=T)
        paths = [fill(one(Complex{T}), 2, 3), fill(one(Complex{T}), 2, 3)]
        coherent = Optics.prepare_spectral_synthesis(paths, ones(T, 1, 1),
            ones(Complex{T}, 1, 2), [ones(Complex{T}, 2, 2)], channel)
        incoherent = Optics.prepare_spectral_synthesis(paths, ones(T, 1, 1),
            ones(Complex{T}, 1, 2), [Matrix{Complex{T}}(I, 2, 2)], channel)
        @test Optics.synthesize_spectrum!(coherent) ≈ fill(T(4), 2, 3)
        @test all(==(T(2)), Optics.synthesize_spectrum!(incoherent))
        fill!(paths[2], -one(Complex{T}))
        @test all(>=(zero(T)), Optics.synthesize_spectrum!(coherent))
        @test maximum(Optics.synthesize_spectrum!(coherent)) <= 16eps(T)
        @test all(==(T(2)), Optics.synthesize_spectrum!(incoherent))

        # The opposite interference phases cancel across two independent
        # wavelength channels. Coherent combination across channels would fail.
        two_channels = SpectralBundle(T[600e-9, 800e-9], T[1, 1]; T=T)
        fill!(paths[2], one(Complex{T}))
        cancellation = Optics.prepare_spectral_synthesis(paths, ones(T, 2, 1),
            Complex{T}[1 1; 1 -1], [ones(Complex{T}, 2, 2) for _ in 1:2],
            two_channels)
        @test Optics.synthesize_spectrum!(cancellation) ≈ fill(T(2), 2, 3)
        @test cancellation.plan.uniform_average

        three_channels = SpectralBundle(T[600e-9, 700e-9, 800e-9], ones(T, 3); T=T)
        uniform_envelope = fill(Complex{T}(0.3, 0.7), 2, 3)
        uniform_coefficients = reshape(Complex{T}[0.15im, 0.8 + im, 1.2im], 3, 1)
        uniform = Optics.prepare_spectral_synthesis([uniform_envelope], ones(T, 3, 1),
            uniform_coefficients, [ones(Complex{T}, 1, 1) for _ in 1:3], three_channels;
            summation=:path_pairs)
        # Uniform averaging retains termwise division by the channel count,
        # including non-power-of-two counts, as the physical15 calibration does.
        divided_intensity = sum(abs2(uniform_envelope[1] * uniform_coefficients[c]) /
            T(3) for c in 1:3)
        @test uniform.plan.uniform_average
        @test all(==(divided_intensity), Optics.synthesize_spectrum!(uniform))

        # A complex off-diagonal coherence and nontrivial interpolation test
        # the sign/conjugation convention against the dense quadratic form.
        anchors = [fill(Complex{T}(k, (-1)^k), 2, 3) for k in 1:6]
        interpolation = T[0.25 0.5 0.25; -0.25 1 0.25]
        gain = Complex{T}[1 im; 2im 0.5]
        coherence = [Complex{T}[1 0.5im; -0.5im 1],
            Complex{T}[1 -0.75; -0.75 1]]
        throughput = T[0.5, 2]
        for summation in (:coherent_modes, :path_pairs)
            mixed = Optics.prepare_spectral_synthesis(anchors, interpolation, gain,
                coherence, two_channels; response=throughput, summation=summation)
            @test Optics.synthesize_spectrum!(mixed) ≈ spectral_synthesis_oracle(
                anchors, interpolation, gain, coherence, two_channels, throughput)
            @test mixed.plan.summation === summation
            @test isconcretetype(eltype(mixed.plan.mode_factors))
            if !coverage_instrumented()
                @test spectral_synthesis_allocations(mixed) == 0
            end
        end

        # An exact PSD matrix and nearly cancelling Float32 fields reproduce
        # the negative intensity from signed pair arithmetic. Squared coherent
        # modes stay nonnegative without clamping the output.
        destructive_values = Complex{T}[0.3 + 0.7im, 0.2 - 0.1im, -0.5 - 0.6im]
        destructive_paths = [fill(value, 2, 3) for value in destructive_values]
        destructive = Optics.prepare_spectral_synthesis(destructive_paths,
            ones(T, 1, 1), ones(Complex{T}, 1, 3), [ones(Complex{T}, 3, 3)], channel)
        @test destructive.plan.summation === :coherent_modes
        @test all(>=(zero(T)), Optics.synthesize_spectrum!(destructive))
        # Near a null, factor/output rounding contributes at squared output
        # precision, while retained positive eigensolve error is at Float64
        # precision. This bound rejects the old Float32 positive-eigen artifact.
        null_tolerance = 64 * (eps(Float64) + Float64(eps(T))^2) *
            sum(abs2, destructive_values)
        @test maximum(destructive.products.output) <= null_tolerance
        if T === Float32
            pairs = Optics.prepare_spectral_synthesis(destructive_paths,
                ones(T, 1, 1), ones(Complex{T}, 1, 3), [ones(Complex{T}, 3, 3)], channel;
                summation=:path_pairs)
            @test all(<(zero(T)), Optics.synthesize_spectrum!(pairs))
        end

        # A tiny negative eigenvalue admitted by the calibration tolerance is
        # projected only in mode factors. The supplied visibility snapshot is
        # retained for reference arithmetic, and final intensity is not clamped.
        near_psd = Complex{T}[1 1 + 2eps(T); 1 + 2eps(T) 1]
        projected = Optics.prepare_spectral_synthesis(paths, ones(T, 1, 1),
            Complex{T}[1 -1], [near_psd], channel)
        factor = Array(projected.plan.mode_factors[1])
        effective = factor * factor'
        @test minimum(eigvals(Hermitian(near_psd))) < zero(T)
        @test projected.plan.visibility[1] == near_psd
        @test effective ≈ near_psd
        @test minimum(eigvals(Hermitian(effective))) >= -8eps(T)
        @test all(>=(zero(T)), Optics.synthesize_spectrum!(projected))
        @test maximum(projected.products.output) <= 16eps(T)
        unprojected = Optics.prepare_spectral_synthesis(paths, ones(T, 1, 1),
            Complex{T}[1 -1], [near_psd], channel; summation=:path_pairs)
        @test all(<(zero(T)), Optics.synthesize_spectrum!(unprojected))
        if !coverage_instrumented()
            @test spectral_synthesis_allocations(destructive) == 0
        end
    end

    T = Float32
    channels = SpectralBundle(T[700e-9], T[1]; T=T)
    envelopes = [ones(Complex{T}, 2, 2), ones(Complex{T}, 2, 2)]
    basis = ones(T, 1, 1)
    coefficients = ones(Complex{T}, 1, 2)
    visibility = [Matrix{Complex{T}}(I, 2, 2)]
    prepare(; e=envelopes, b=basis, c=coefficients, v=visibility, kwargs...) =
        Optics.prepare_spectral_synthesis(e, b, c, v, channels; kwargs...)
    @test_throws DimensionMismatchError prepare(b=ones(T, 2, 1))
    @test_throws DimensionMismatchError prepare(e=envelopes[1:1])
    @test_throws DimensionMismatchError prepare(output=zeros(T, 3, 2))
    @test_throws DimensionMismatchError prepare(v=[ones(Complex{T}, 1, 1)])
    @test_throws InvalidConfiguration prepare(e=[envelopes[1], envelopes[1]])
    aliased_output = reshape(reinterpret(T, vec(envelopes[1])), 2, 4)
    aliased_envelope = [reshape(envelopes[1], 1, 4)]
    @test_throws InvalidConfiguration Optics.prepare_spectral_synthesis(
        aliased_envelope, basis, ones(Complex{T}, 1, 1),
        [ones(Complex{T}, 1, 1)], channels;
        output=view(aliased_output, 1:1, 1:4))
    @test_throws InvalidConfiguration prepare(response=T[-1])
    @test_throws InvalidConfiguration prepare(response=[-floatmin(Float64)])
    @test_throws InvalidConfiguration prepare(response=T[Inf])
    @test_throws InvalidConfiguration prepare(response=one(T))
    @test_throws InvalidConfiguration prepare(summation=:unknown)
    @test_throws InvalidConfiguration prepare(summation=missing)
    @test_throws InvalidConfiguration prepare(b=fill(Float64(typemax(T)) * 2, 1, 1))
    @test_throws InvalidConfiguration prepare(c=fill(ComplexF64(Inf), 1, 2))
    @test_throws InvalidConfiguration prepare(v=[Complex{T}[1 1; 0 1]])
    @test_throws InvalidConfiguration prepare(v=[Complex{T}[2 0; 0 1]])
    @test_throws InvalidConfiguration prepare(v=[Complex{T}[1 2; 2 1]])
    @test_throws InvalidConfiguration Optics.prepare_spectral_synthesis(
        [copy(envelopes[1]) for _ in 1:3], basis, ones(Complex{T}, 1, 3),
        [Complex{T}[1 0.9 0.9; 0.9 1 -0.9; 0.9 -0.9 1]], channels)
    @test_throws InvalidConfiguration prepare(v=[Complex{T}[1 NaN; NaN 1]])
    @test_throws InvalidConfiguration prepare(b=ones(Complex{T}, 1, 1))
    @test_throws InvalidConfiguration prepare(c=ones(T, 1, 2))
    @test_throws InvalidConfiguration prepare(output=zeros(Float64, 2, 2))
    @test_throws InvalidConfiguration prepare(output=zeros(T, 4))
    foreign = SpectralSynthesisForeignMatrix.(envelopes)
    @test_throws InvalidConfiguration prepare(e=foreign, output=zeros(T, 2, 2))
    shifted = SpectralSynthesisShiftedMatrix.(envelopes)
    @test_throws ArgumentError prepare(e=shifted, output=zeros(T, 2, 2))
    @test_throws InvalidConfiguration Optics.prepare_spectral_synthesis(
        Matrix{Complex{T}}[], basis, coefficients, visibility, channels)
    # Nonoverlapping views are admitted when mightalias proves disjoint storage.
    parent = ones(Complex{T}, 4, 2)
    viewed = prepare(e=[view(parent, 1:2, :), view(parent, 3:4, :)])
    @test all(==(T(2)), Optics.synthesize_spectrum!(viewed))
    @test_throws InvalidConfiguration prepare(e=[view(parent, 1:2, :), view(parent, 2:3, :)])

    for name in (:SpectralSynthesisPlan, :SpectralSynthesisWorkspace,
            :SpectralSynthesisProducts, :PreparedSpectralSynthesis,
            :prepare_spectral_synthesis, :synthesize_spectrum!)
        @test Base.ispublic(Optics, name)
        @test !Base.isexported(Optics, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
end

"""Run-immutable channel quadrature, interpolation, and path-coherence calibration."""
struct SpectralSynthesisPlan{T,B,C,V,M,S,R,W,D,E,K}
    shape::Tuple{Int,Int}
    n_anchor::Int
    n_path::Int
    basis::B
    coefficients::C
    visibility::V
    mode_factors::M
    summation::Symbol
    channels::S
    response::R
    weights::W
    uniform_average::Bool
    channel_count::T
    backend::E
    device::D
    resident_calibration::K
end

"""Replaceable interpolated field and accelerator packing scratch for spectral synthesis."""
struct SpectralSynthesisWorkspace{F,A,P,C}
    channel_fields::F
    coherent_field::A
    packed_envelopes::P
    packed_channel_fields::C
end

# Resident copies are immutable calibration; packed detector fields remain scratch.
struct _SpectralSynthesisResidentCalibration{B,C,V,M,W}
    basis::B
    coefficients::C
    visibility::V
    mode_factors::M
    weights::W
end

_spectral_synthesis_resident_calibration(::ScalarCPUStyle, output, basis,
    coefficients, visibility, mode_factors, weights) = nothing

function _spectral_synthesis_resident_calibration(::AcceleratorStyle, output,
        basis, coefficients, visibility, mode_factors, weights)
    function resident(array)
        # FixedSizeArray concatenation retains fixed-size host storage, which
        # accelerator copy methods do not necessarily upload without scalar access.
        host = Array(array)
        return copyto!(similar(output, eltype(host), size(host)), host)
    end
    return _SpectralSynthesisResidentCalibration(resident(basis),
        resident(coefficients), resident(cat(visibility...; dims=3)),
        resident(cat(mode_factors...; dims=3)), resident(weights))
end

_spectral_synthesis_packed_scratch(::ScalarCPUStyle, envelope, n_path,
    n_anchor, n_channel) = (nothing, nothing)

function _spectral_synthesis_packed_scratch(::AcceleratorStyle, envelope,
        n_path, n_anchor, n_channel)
    return (similar(envelope, length(envelope), n_anchor * n_path),
        similar(envelope, length(envelope), n_path, n_channel))
end

"""Caller-visible intensity output of spectral synthesis."""
struct SpectralSynthesisProducts{O}
    output::O
end

struct _SpectralSynthesisPreparationToken end
const _SPECTRAL_SYNTHESIS_PREPARATION_TOKEN = _SpectralSynthesisPreparationToken()

"""
    PreparedSpectralSynthesis

Single-writer owner of exact path-envelope bindings, frozen spectral calibration,
replaceable field scratch, and a caller-visible intensity output. Construct with
`prepare_spectral_synthesis`; input values may change between completed calls.
"""
struct PreparedSpectralSynthesis{I,P,W,O}
    envelopes::I
    plan::P
    workspace::W
    products::O

    function PreparedSpectralSynthesis(::_SpectralSynthesisPreparationToken,
            envelopes::I, plan::P, workspace::W, products::O) where {I,P,W,O}
        return new{I,P,W,O}(envelopes, plan, workspace, products)
    end
end

function _spectral_synthesis_default_output(envelopes)
    isempty(envelopes) && throw(InvalidConfiguration(
        "spectral synthesis requires at least one path envelope"))
    envelope = first(envelopes)
    return _with_compute_device(compute_device(envelope)) do
        similar(envelope, real(eltype(envelope)))
    end
end

function _spectral_synthesis_calibration(array, ::Type{T}, label) where {T}
    host = T.(Array(array))
    all(isfinite, host) || throw(InvalidConfiguration(
        "$label must remain finite at the output precision"))
    return host
end

function _spectral_synthesis_visibility(matrix, ::Type{T}, n_path) where {T}
    eltype(matrix) <: Complex || throw(InvalidConfiguration(
        "spectral synthesis visibility must have complex elements"))
    Base.require_one_based_indexing(matrix)
    size(matrix) == (n_path, n_path) || throw(DimensionMismatchError(
        "spectral synthesis visibility must have one row and column per path"))
    host = _spectral_synthesis_calibration(matrix, Complex{T},
        "spectral synthesis visibility")
    tolerance = 8eps(T)
    for k in 1:n_path
        abs(host[k, k] - one(T)) <= tolerance || throw(InvalidConfiguration(
            "spectral synthesis visibility must have unit diagonal"))
        for l in (k + 1):n_path
            abs(host[k, l] - conj(host[l, k])) <= tolerance ||
                throw(InvalidConfiguration(
                    "spectral synthesis visibility must be Hermitian"))
        end
    end
    # Cold Float64 factorization avoids output-precision positive eigenvalue
    # error becoming spurious power in a rank-deficient Float32 coherence.
    decomposition = eigen(Hermitian(ComplexF64.(host)))
    all(isfinite, decomposition.values) &&
        minimum(decomposition.values) >= -tolerance * n_path ||
        throw(InvalidConfiguration(
            "spectral synthesis visibility must be positive semidefinite"))
    # Only eigenvalues already admitted within the PSD roundoff tolerance are
    # projected to zero. The resulting factor represents V ≈ F * F'.
    factor = decomposition.vectors *
        Diagonal(sqrt.(max.(decomposition.values, zero(T))))
    return (visibility=FixedSizeMatrixDefault{Complex{T}}(host),
        mode_factor=FixedSizeMatrixDefault{Complex{T}}(factor))
end

"""
    prepare_spectral_synthesis(envelopes, basis, coefficients, visibility, channels;
                              response=ones(length(channels)), output=...,
                              summation=:coherent_modes)

Prepare wavelength-channel synthesis on one caller-declared physical detector
grid. `envelopes` holds `ComplexF32` or `ComplexF64` matrices in path-major,
anchor-minor order. All bound matrices have equal nonempty shape, one-based axes,
one exact compute device and backend, and no aliases. Their concrete storage
types must match. The output has the corresponding real precision; its default
is newly allocated storage similar to the first envelope.

The real interpolation `basis` has shape `(n_channel, n_anchor)` and the complex
path `coefficients` have shape `(n_channel, n_path)`. Each complex
`visibility[channel]` is a unit-diagonal Hermitian positive semidefinite path
coherence matrix of shape `(n_path, n_path)`. Hermitian and diagonal tolerances
are eight machine epsilons; the minimum eigenvalue tolerance is that value times
the path count. Calibration is copied into fixed-size host storage, with resident
accelerator copies when needed, and must remain finite at output precision. It
must not alias any bound matrix.

`summation=:coherent_modes` is the default. Preparation factors each coherence
matrix in Float64 as `V ≈ F * F'` and freezes factors at output precision,
projecting only eigenvalues admitted as negative roundoff
within the PSD tolerance to zero. Execution forms each coherent mode as
`sum(field[k] * F[k, mode])` and adds its squared magnitude. Contributions are
nonnegative by construction; the effective coherence matrix can differ from
the supplied matrix by the tolerated eigenvalue projection and factorization
roundoff. `summation=:path_pairs` preserves the supplied diagonal/pair arithmetic
and termwise division for uniform channel weights and unit response. That
reference policy can produce small negative intensities at destructive
interference because signed terms cancel at finite precision. Neither policy
clamps the final intensity.

`channels` is a `SpectralBundle` whose normalized weights are photon-number
fractions. `response` is finite nonnegative multiplicative throughput or quantum
efficiency per channel, and is not renormalized. Each channel combines path
fields according to its coherence before incoherent addition across channels.
The result has
units of squared input amplitude times channel weights and response. The caller
supplies calibrated coefficients and envelopes already on the same physical
grid; this operation performs no propagation, resampling, or inferred photometry.

Preparation may allocate and observe calibration on the host. Execution queues
backend kernels on the caller's current accelerator stream, with no completion
boundary; the caller owns completion before reading output or changing inputs.
Accelerators pack the live envelopes each call, interpolate all path/channel
fields, then accumulate channels and coherent modes or path pairs in their
original order. Packed scratch adds `length(output) * n_path *
(n_anchor + length(channels))` complex elements. CPU execution keeps the original
per-path scratch and independent broadcast arithmetic.
"""
function prepare_spectral_synthesis(envelopes::AbstractVector{<:AbstractMatrix},
        basis::AbstractMatrix, coefficients::AbstractMatrix,
        visibility::AbstractVector{<:AbstractMatrix}, channels::SpectralBundle;
        response=ones(length(channels)),
        output=_spectral_synthesis_default_output(envelopes),
        summation=:coherent_modes)
    (summation === :coherent_modes || summation === :path_pairs) || throw(InvalidConfiguration(
        "spectral synthesis summation must be :coherent_modes or :path_pairs"))
    isempty(envelopes) && throw(InvalidConfiguration(
        "spectral synthesis requires at least one path envelope"))
    output isa AbstractMatrix || throw(InvalidConfiguration(
        "spectral synthesis output must be a matrix"))
    response isa AbstractVector && eltype(response) <: Real ||
        throw(InvalidConfiguration("spectral synthesis response must be a real vector"))
    Base.require_one_based_indexing(envelopes, basis, coefficients, visibility,
        channels.samples, response, output)
    T = eltype(output)
    T in (Float32, Float64) || throw(InvalidConfiguration(
        "spectral synthesis output must have Float32 or Float64 elements"))
    eltype(basis) <: Real || throw(InvalidConfiguration(
        "spectral synthesis basis must have real elements"))
    eltype(coefficients) <: Complex || throw(InvalidConfiguration(
        "spectral synthesis coefficients must have complex elements"))
    n_channel = length(channels)
    n_anchor = size(basis, 2)
    n_path = size(coefficients, 2)
    n_anchor > 0 && n_path > 0 || throw(InvalidConfiguration(
        "spectral synthesis requires at least one anchor and path"))
    n_anchor <= typemax(Int) ÷ n_path || throw(InvalidConfiguration(
        "spectral synthesis envelope count overflows Int"))
    size(basis, 1) == n_channel && size(coefficients, 1) == n_channel &&
        length(visibility) == n_channel && length(response) == n_channel ||
        throw(DimensionMismatchError("spectral synthesis calibration channel counts differ"))
    length(envelopes) == n_anchor * n_path || throw(DimensionMismatchError(
        "spectral synthesis envelope count must equal anchors times paths"))
    shape = size(output)
    all(>(0), shape) || throw(InvalidConfiguration(
        "spectral synthesis detector grid must be nonempty"))
    device = compute_device(output)
    envelope_type = typeof(first(envelopes))
    isconcretetype(envelope_type) || throw(InvalidConfiguration(
        "spectral synthesis envelopes must have concrete storage"))
    for (i, envelope) in pairs(envelopes)
        Base.require_one_based_indexing(envelope)
        typeof(envelope) === envelope_type && eltype(envelope) === Complex{T} ||
            throw(InvalidConfiguration(
                "spectral synthesis envelopes must share concrete storage and output precision"))
        size(envelope) == shape || throw(DimensionMismatchError(
            "spectral synthesis envelopes must share the output shape"))
        require_same_backend(output, envelope)
        compute_device(envelope) == device || _throw_compute_device_error(
            :prepare_spectral_synthesis, :wrong_device, device,
            "spectral synthesis envelope occupies $(compute_device(envelope))")
        Base.mightalias(envelope, output) && throw(InvalidConfiguration(
            "spectral synthesis envelopes and output must not alias"))
        for j in 1:(i - 1)
            Base.mightalias(envelope, envelopes[j]) && throw(InvalidConfiguration(
                "spectral synthesis path envelopes must not alias"))
        end
    end
    for calibration in (basis, coefficients, response, visibility...)
        (Base.mightalias(output, calibration) ||
            any(envelope -> Base.mightalias(envelope, calibration), envelopes)) &&
            throw(InvalidConfiguration(
                "spectral synthesis calibration must not alias bound storage"))
    end
    _validated_spectral_weight_sum(channels.samples)
    host_basis = _spectral_synthesis_calibration(basis, T, "spectral synthesis basis")
    host_coefficients = _spectral_synthesis_calibration(coefficients, Complex{T},
        "spectral synthesis coefficients")
    supplied_response = Array(response)
    all(value -> isfinite(value) && value >= zero(value), supplied_response) ||
        throw(InvalidConfiguration(
            "spectral synthesis response must be finite and nonnegative"))
    host_response = _spectral_synthesis_calibration(supplied_response, T,
        "spectral synthesis response")
    host_channels = [SpectralSample{T}(T(sample.wavelength), T(sample.weight))
        for sample in channels]
    _validated_spectral_weight_sum(host_channels)
    isapprox(sum(sample.weight for sample in host_channels), one(T);
        atol=sqrt(eps(T)), rtol=sqrt(eps(T))) || throw(InvalidConfiguration(
            "spectral synthesis photon-number weights must sum to one"))
    weights = [host_channels[c].weight * host_response[c] for c in 1:n_channel]
    all(isfinite, weights) || throw(InvalidConfiguration(
        "spectral synthesis weighted response overflows output precision"))
    coherence = [_spectral_synthesis_visibility(v, T, n_path) for v in visibility]
    host_visibility = [value.visibility for value in coherence]
    mode_factors = [value.mode_factor for value in coherence]
    uniform_average = all(sample -> sample.weight == inv(T(n_channel)), host_channels) &&
        all(==(one(T)), host_response)
    style = execution_style(output)
    resident_calibration, packed_envelopes, packed_channel_fields =
        _with_compute_device(device) do
            calibration = _spectral_synthesis_resident_calibration(style, output,
                host_basis, host_coefficients, host_visibility, mode_factors, weights)
            packed = _spectral_synthesis_packed_scratch(style, first(envelopes),
                n_path, n_anchor, n_channel)
            return calibration, packed...
        end
    plan = SpectralSynthesisPlan(shape, n_anchor, n_path,
        FixedSizeMatrixDefault{T}(host_basis),
        FixedSizeMatrixDefault{Complex{T}}(host_coefficients),
        FixedSizeVectorDefault{eltype(host_visibility)}(host_visibility),
        FixedSizeVectorDefault{eltype(mode_factors)}(mode_factors), summation,
        FixedSizeVectorDefault{SpectralSample{T}}(host_channels),
        FixedSizeVectorDefault{T}(host_response), FixedSizeVectorDefault{T}(weights),
        uniform_average, T(n_channel), backend(output), device, resident_calibration)
    fields, coherent_field = _with_compute_device(device) do
        return [similar(first(envelopes)) for _ in 1:n_path], similar(first(envelopes))
    end
    return PreparedSpectralSynthesis(_SPECTRAL_SYNTHESIS_PREPARATION_TOKEN,
        FixedSizeVectorDefault{envelope_type}(envelopes), plan,
        SpectralSynthesisWorkspace(FixedSizeVectorDefault{eltype(fields)}(fields),
            coherent_field, packed_envelopes, packed_channel_fields),
        SpectralSynthesisProducts(output))
end

"""
    synthesize_spectrum!(owner::PreparedSpectralSynthesis)

Synthesize the bound path envelopes into `owner.products.output`, returning that
exact matrix. Accelerator work is submitted on the caller's current stream on
the prepared device without synchronization or host observation. Warmed CPU
execution allocates no heap storage. All bindings and calibration stay fixed.
"""
function synthesize_spectrum!(owner::PreparedSpectralSynthesis)
    return _synthesize_spectrum!(execution_style(owner.products.output), owner)
end

_synthesize_spectrum!(::ScalarCPUStyle, owner) = _synthesize_spectrum_reference!(owner)

# Keep the original broadcast executor callable on either backend as an independent
# arithmetic and performance reference. It uses only the original per-path scratch.
function _synthesize_spectrum_reference!(owner::PreparedSpectralSynthesis)
    plan = owner.plan
    fields = owner.workspace.channel_fields
    out = owner.products.output
    fill!(out, zero(eltype(out)))
    for c in eachindex(plan.channels)
        for k in 1:plan.n_path
            field = fields[k]
            fill!(field, zero(eltype(field)))
            for a in 1:plan.n_anchor
                basis = plan.basis[c, a]
                envelope = owner.envelopes[(k - 1) * plan.n_anchor + a]
                @. field += basis * envelope
            end
            coefficient = plan.coefficients[c, k]
            field .*= coefficient
        end
        weight = plan.weights[c]
        divisor = plan.channel_count
        if plan.summation === :coherent_modes
            coherent_field = owner.workspace.coherent_field
            factor = plan.mode_factors[c]
            for mode in 1:plan.n_path
                fill!(coherent_field, zero(eltype(coherent_field)))
                for k in 1:plan.n_path
                    field = fields[k]
                    coefficient = factor[k, mode]
                    @. coherent_field += field * coefficient
                end
                if plan.uniform_average
                    @. out += abs2(coherent_field) / divisor
                else
                    @. out += abs2(coherent_field) * weight
                end
            end
            continue
        end
        for k in 1:plan.n_path
            field = fields[k]
            diagonal = real(plan.visibility[c][k, k])
            if plan.uniform_average
                @. out += diagonal * abs2(field) / divisor
            else
                @. out += diagonal * abs2(field) * weight
            end
            for l in (k + 1):plan.n_path
                other = fields[l]
                cross = plan.visibility[c][k, l]
                if plan.uniform_average
                    @. out += 2 * real(field * cross * conj(other)) / divisor
                else
                    @. out += 2 * real(field * cross * conj(other)) * weight
                end
            end
        end
    end
    return out
end

# Packing observes the exact bound matrices on every call, including strided views.
@kernel function _spectral_synthesis_pack_kernel!(packed, envelope, envelope_index)
    row, column = @index(Global, NTuple)
    if row <= size(envelope, 1) && column <= size(envelope, 2)
        pixel = row + (column - 1) * size(envelope, 1)
        packed[pixel, envelope_index] = envelope[row, column]
    end
end

@kernel function _spectral_synthesis_fields_kernel!(fields, envelopes, basis,
        coefficients, n_anchor)
    pixel, path, channel = @index(Global, NTuple)
    if pixel <= size(fields, 1) && path <= size(fields, 2) && channel <= size(fields, 3)
        field = zero(eltype(fields))
        for anchor in 1:n_anchor
            field += basis[channel, anchor] *
                envelopes[pixel, (path - 1) * n_anchor + anchor]
        end
        fields[pixel, path, channel] = field * coefficients[channel, path]
    end
end

@inline function _spectral_synthesis_channel(intensity, field_power, weight,
        divisor, uniform_average)
    if uniform_average
        return intensity + field_power / divisor
    end
    return intensity + field_power * weight
end

@inline function _spectral_synthesis_accumulate(intensity, fields, factors,
        pixel, channel, weight, divisor, uniform_average, ::Val{:coherent_modes})
    for mode in axes(fields, 2)
        coherent = zero(eltype(fields))
        for path in axes(fields, 2)
            coherent += fields[pixel, path, channel] * factors[path, mode, channel]
        end
        intensity = _spectral_synthesis_channel(intensity, abs2(coherent),
            weight, divisor, uniform_average)
    end
    return intensity
end

@inline function _spectral_synthesis_accumulate(intensity, fields, visibility,
        pixel, channel, weight, divisor, uniform_average, ::Val{:path_pairs})
    for path in axes(fields, 2)
        field = fields[pixel, path, channel]
        diagonal = real(visibility[path, path, channel])
        intensity = _spectral_synthesis_channel(intensity, diagonal * abs2(field),
            weight, divisor, uniform_average)
        for other_path in (path + 1):size(fields, 2)
            other = fields[pixel, other_path, channel]
            cross = visibility[path, other_path, channel]
            intensity = _spectral_synthesis_channel(intensity,
                2 * real(field * cross * conj(other)), weight, divisor,
                uniform_average)
        end
    end
    return intensity
end

@kernel function _spectral_synthesis_accumulate_kernel!(output, fields,
        coherence, weights, divisor, uniform_average, summation)
    pixel = @index(Global, Linear)
    if pixel <= length(output)
        intensity = zero(eltype(output))
        for channel in axes(fields, 3)
            intensity = _spectral_synthesis_accumulate(intensity, fields, coherence,
                pixel, channel, weights[channel], divisor, uniform_average, summation)
        end
        output[pixel] = intensity
    end
end

function _synthesize_spectrum!(style::AcceleratorStyle, owner)
    return _synthesize_spectrum!(style, owner, owner.plan.resident_calibration,
        owner.workspace.packed_envelopes, owner.workspace.packed_channel_fields)
end

function _synthesize_spectrum!(style::AcceleratorStyle, owner, calibration,
        packed_envelopes, fields)
    plan = owner.plan
    for envelope_index in eachindex(owner.envelopes)
        launch_kernel_async!(style, _spectral_synthesis_pack_kernel!,
            packed_envelopes, owner.envelopes[envelope_index], envelope_index;
            ndrange=plan.shape)
    end
    launch_kernel_async!(style, _spectral_synthesis_fields_kernel!, fields,
        packed_envelopes, calibration.basis, calibration.coefficients,
        plan.n_anchor; ndrange=size(fields))
    if plan.summation === :coherent_modes
        _spectral_synthesis_accumulate!(style, owner, fields, calibration.mode_factors,
            calibration.weights, Val(:coherent_modes))
    else
        _spectral_synthesis_accumulate!(style, owner, fields, calibration.visibility,
            calibration.weights, Val(:path_pairs))
    end
    return owner.products.output
end

function _spectral_synthesis_accumulate!(style, owner, fields, coherence, weights,
        summation)
    plan = owner.plan
    launch_kernel_async!(style, _spectral_synthesis_accumulate_kernel!,
        owner.products.output, fields, coherence, weights,
        plan.channel_count, plan.uniform_average,
        summation; ndrange=length(owner.products.output))
    return nothing
end

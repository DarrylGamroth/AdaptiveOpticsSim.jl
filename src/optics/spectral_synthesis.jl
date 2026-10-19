"""Run-immutable channel quadrature, interpolation, and path-coherence calibration."""
struct SpectralSynthesisPlan{T,B,C,V,M,S,R,W,D,E}
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
end

"""Replaceable per-path interpolated field scratch for spectral synthesis."""
struct SpectralSynthesisWorkspace{F,A}
    channel_fields::F
    coherent_field::A
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
the path count. Calibration is copied into fixed-size host storage and must
remain finite at output precision. It must not alias any bound matrix.

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
broadcast work on the caller's current accelerator stream, with no completion
boundary; the caller owns completion before reading output or changing inputs.
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
    plan = SpectralSynthesisPlan(shape, n_anchor, n_path,
        FixedSizeMatrixDefault{T}(host_basis),
        FixedSizeMatrixDefault{Complex{T}}(host_coefficients),
        FixedSizeVectorDefault{eltype(host_visibility)}(host_visibility),
        FixedSizeVectorDefault{eltype(mode_factors)}(mode_factors), summation,
        FixedSizeVectorDefault{SpectralSample{T}}(host_channels),
        FixedSizeVectorDefault{T}(host_response), FixedSizeVectorDefault{T}(weights),
        uniform_average, T(n_channel), backend(output), device)
    fields, coherent_field = _with_compute_device(device) do
        return [similar(first(envelopes)) for _ in 1:n_path], similar(first(envelopes))
    end
    return PreparedSpectralSynthesis(_SPECTRAL_SYNTHESIS_PREPARATION_TOKEN,
        FixedSizeVectorDefault{envelope_type}(envelopes), plan,
        SpectralSynthesisWorkspace(FixedSizeVectorDefault{eltype(fields)}(fields),
            coherent_field),
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

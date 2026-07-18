#
# Prepared focal-plane modulation
#
# Pyramid and BioEdge sensors may use the same pupil-plane tip/tilt path while
# retaining independent focal-plane masks and propagation workspaces. The
# weights below are optical quadrature weights: they average instantaneous
# intensities over one modulation cycle and never represent elapsed exposure.
#

abstract type AbstractFocalPlaneModulation end

"""A stationary focal-plane sensor with one unit-weight optical sample."""
struct NoModulation <: AbstractFocalPlaneModulation end

"""
    CircularModulation(radius; samples, phase_offset=0)

Uniform circular modulation in units of `lambda / D`. `samples` is the number
of equally weighted points over one cycle and `phase_offset` is in radians.
"""
struct CircularModulation{T<:AbstractFloat} <: AbstractFocalPlaneModulation
    radius::T
    samples::Int
    phase_offset::T

    function CircularModulation(radius::T, samples::Int,
        phase_offset::T) where {T<:AbstractFloat}
        isfinite(radius) && radius >= zero(T) || throw(InvalidConfiguration(
            "circular modulation radius must be finite and nonnegative"))
        samples >= 1 || throw(InvalidConfiguration(
            "circular modulation samples must be >= 1"))
        isfinite(phase_offset) || throw(InvalidConfiguration(
            "circular modulation phase offset must be finite"))
        return new{T}(radius, samples, phase_offset)
    end
end

function CircularModulation(radius::Real; samples::Int,
    phase_offset::Real=0, T::Type{<:AbstractFloat}=Float64)
    return CircularModulation(T(radius), samples, T(phase_offset))
end

"""
    SampledModulation(points; weights=nothing)

An arbitrary focal-plane modulation path. Each point is an `(x, y)` offset in
units of `lambda / D`. Weights are copied, normalized, finite, nonnegative
optical quadrature weights; at least one must be positive.
"""
struct SampledModulation{T<:AbstractFloat} <: AbstractFocalPlaneModulation
    points::Vector{NTuple{2,T}}
    weights::Vector{T}
end

function SampledModulation(points; weights=nothing,
    T::Type{<:AbstractFloat}=Float64)
    isempty(points) && throw(InvalidConfiguration(
        "sampled modulation must contain at least one point"))
    owned_points = Vector{NTuple{2,T}}(undef, length(points))
    @inbounds for index in eachindex(points)
        point = points[index]
        length(point) >= 2 || throw(InvalidConfiguration(
            "each sampled modulation point must contain x and y offsets"))
        x = T(point[1])
        y = T(point[2])
        isfinite(x) && isfinite(y) || throw(InvalidConfiguration(
            "sampled modulation offsets must be finite"))
        owned_points[index] = (x, y)
    end
    raw_weights = if weights === nothing
        inferred = Vector{T}(undef, length(points))
        @inbounds for index in eachindex(points)
            point = points[index]
            inferred[index] = length(point) >= 3 ? T(point[3]) : one(T)
        end
        inferred
    else
        length(weights) == length(points) || throw(InvalidConfiguration(
            "sampled modulation weights must match the point count"))
        T.(collect(weights))
    end
    weight_sum = zero(T)
    @inbounds for weight in raw_weights
        isfinite(weight) && weight >= zero(T) || throw(InvalidConfiguration(
            "sampled modulation weights must be finite and nonnegative"))
        weight_sum += weight
    end
    weight_sum > zero(T) || throw(InvalidConfiguration(
        "sampled modulation requires at least one positive weight"))
    raw_weights ./= weight_sum
    return SampledModulation{T}(owned_points, raw_weights)
end

"""Backend-bound modulation phases and host-resident amplitude weights."""
struct PreparedFocalPlaneModulation{P<:AbstractFocalPlaneModulation,A,V}
    policy::P
    phases::A
    amplitude_weights::V
end

@inline modulation_point_count(::NoModulation) = 1
@inline modulation_point_count(policy::CircularModulation) = policy.samples
@inline modulation_point_count(policy::SampledModulation) = length(policy.points)
@inline modulation_point_count(prepared::PreparedFocalPlaneModulation) =
    size(prepared.phases, 3)

@inline modulation_offset(::NoModulation, ::Int, ::Type{T}) where {T} =
    (zero(T), zero(T))

@inline function modulation_offset(policy::CircularModulation, index::Int,
    ::Type{T}) where {T}
    angle = T(policy.phase_offset) + T(2pi * (index - 1) / policy.samples)
    sine, cosine = sincos(angle)
    radius = T(policy.radius)
    return radius * cosine, radius * sine
end

@inline function modulation_offset(policy::SampledModulation, index::Int,
    ::Type{T}) where {T}
    point = @inbounds policy.points[index]
    return T(point[1]), T(point[2])
end

@inline modulation_weight(::NoModulation, ::Int, ::Type{T}) where {T} = one(T)
@inline modulation_weight(policy::CircularModulation, ::Int,
    ::Type{T}) where {T} = inv(T(policy.samples))
@inline modulation_weight(policy::SampledModulation, index::Int,
    ::Type{T}) where {T} = T(@inbounds policy.weights[index])

function prepare_focal_plane_modulation(policy::AbstractFocalPlaneModulation,
    resolution::Int, backend_storage::AbstractArray, ::Type{T}) where {
    T<:AbstractFloat,
}
    resolution >= 1 || throw(InvalidConfiguration(
        "modulation pupil resolution must be positive"))
    point_count = modulation_point_count(policy)
    phases = similar(backend_storage, Complex{T}, resolution, resolution,
        point_count)
    host_phases = Array{Complex{T}}(undef, resolution, resolution, point_count)
    amplitude_weights = Vector{T}(undef, point_count)
    coordinates = range(-T(pi), T(pi); length=resolution)
    @inbounds for point_index in 1:point_count
        offset_x, offset_y = modulation_offset(policy, point_index, T)
        amplitude_weights[point_index] = sqrt(modulation_weight(policy,
            point_index, T))
        for axis_2 in 1:resolution, axis_1 in 1:resolution
            phase = offset_x * coordinates[axis_2] +
                offset_y * coordinates[axis_1]
            host_phases[axis_1, axis_2, point_index] = cis(phase)
        end
    end
    copyto!(phases, host_phases)
    return PreparedFocalPlaneModulation(policy, phases, amplitude_weights)
end

function legacy_modulation_policy(modulation::T,
    modulation_points::Union{Int,Nothing}, extra_modulation_factor::Int,
    delta_theta::T, user_modulation_path) where {T<:AbstractFloat}
    if user_modulation_path !== nothing
        return SampledModulation(user_modulation_path; T=T)
    end
    if iszero(modulation)
        return NoModulation()
    end
    samples = if modulation_points === nothing
        perimeter = T(2pi) * modulation
        max(1, 4 * Int(extra_modulation_factor + ceil(perimeter / 4)))
    else
        modulation_points >= 1 || throw(InvalidConfiguration(
            "modulation_points must be >= 1"))
        modulation_points
    end
    return CircularModulation(modulation, samples, delta_theta)
end

@inline calibration_modulation_policy(policy::SampledModulation,
    radius, phase_offset) = policy

function calibration_modulation_policy(
    policy::Union{NoModulation,CircularModulation}, radius::T,
    phase_offset::T) where {T<:AbstractFloat}
    iszero(radius) && return NoModulation()
    return CircularModulation(radius, modulation_point_count(policy),
        phase_offset)
end

@inline _modulated_input_storage(input::PupilFunction) = input.opd
@inline _modulated_input_storage(input::ElectricField) = input.values

function require_modulated_wfs_input(input::PupilFunction)
    validate_wfs_optical_input(input)
    input.metadata.coordinate_domain isa MetricCoordinates ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "modulated WFS pupil input must use metric coordinates"))
    input.metadata.normalization isa DimensionlessNormalization ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "PupilFunction amplitude must be dimensionless"))
    input.metadata.spatial_measure isa PointSampledMeasure ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "PupilFunction amplitude must be point sampled"))
    input.metadata.coherence isa CoherentFieldCombination ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "modulated WFS pupil input must be coherent"))
    input.metadata.spectral isa AchromaticSpectralCoordinate ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "PupilFunction input must be achromatic"))
    return input
end

function require_modulated_wfs_input(input::ElectricField)
    validate_wfs_optical_input(input)
    input.metadata.coordinate_domain isa MetricCoordinates ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "modulated WFS electric-field input must use metric coordinates"))
    input.metadata.normalization isa PhotonRateNormalization ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "electric-field input must carry photon-rate normalization"))
    input.metadata.spatial_measure isa CellIntegratedMeasure ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "electric-field input must carry cell-integrated photon rate"))
    input.metadata.coherence isa CoherentFieldCombination ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "modulated WFS electric-field input must be coherent"))
    input.metadata.spectral isa MonochromaticChannel ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "electric-field input must declare one monochromatic channel"))
    return input
end

@inline function modulated_input_wavelength(input::ElectricField)
    return input.metadata.spectral.wavelength_m
end

function modulated_input_wavelength(input::PupilFunction, source)
    source === nothing && throw(WFSPreparationError(:optical_formation,
        :radiometry, "PupilFunction formation requires an illumination source"))
    return wavelength(source)
end

@inline modulated_wfs_propagation_storage(front_end) =
    front_end.propagation.field

function require_modulated_wfs_domains(front_end, input, output::IntensityMap)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "modulated WFS input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:optical_formation, :device,
            "modulated WFS input and output occupy different devices"))
    storage = modulated_wfs_propagation_storage(front_end)
    typeof(input.metadata.backend) === typeof(backend(storage)) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "modulated WFS input and propagation backends differ"))
    input.metadata.device == plane_device(storage) ||
        throw(WFSPreparationError(:optical_formation, :device,
            "modulated WFS input and propagation occupy different devices"))
    return nothing
end

function require_four_pupil_rate_map(output::IntensityMap, expected_size,
    wavelength_m)
    validate_wfs_optical_products(output)
    output.metadata.coordinate_domain isa NormalizedPupilCoordinates ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "four-pupil detector output must use normalized pupil coordinates"))
    output.metadata.spatial_measure isa CellIntegratedMeasure ||
        throw(WFSPreparationError(:optical_formation, :radiometry,
            "four-pupil detector output must carry cell-integrated rate"))
    size(output.values) == expected_size || throw(WFSPreparationError(
        :optical_formation, :shape,
        "four-pupil detector output has the wrong prepared dimensions"))
    channel = output.metadata.spectral
    channel isa MonochromaticChannel &&
        channel.wavelength_m == wavelength_m ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "four-pupil detector output wavelength does not match its input"))
    return output
end

struct FourPupilSpectralComponent{S,T<:AbstractFloat} <: AbstractSource
    source::S
    wavelength_m::T
    photon_rate_m2_s::T
end

@inline wavelength(source::FourPupilSpectralComponent) = source.wavelength_m
@inline photon_irradiance(source::FourPupilSpectralComponent) =
    source.photon_rate_m2_s

# LGS image formation is prepared at the optical-stage boundary. These values
# own only immutable execution data; the front-end propagation workspace
# remains the single writer of FFT and scratch storage.
abstract type AbstractPreparedFourPupilLGS end

struct NoPreparedFourPupilLGS <: AbstractPreparedFourPupilLGS end

struct PreparedFourPupilElongation{K<:AbstractVector} <:
    AbstractPreparedFourPupilLGS
    kernel::K
    half_width::Int
end

struct PreparedFourPupilSodiumProfile{K<:AbstractMatrix} <:
    AbstractPreparedFourPupilLGS
    kernel_fft::K
end

@inline _four_pupil_lgs_source(source) = nothing
@inline _four_pupil_lgs_source(source::LGSSource) = source
@inline _four_pupil_lgs_source(
    source::FourPupilSpectralComponent{<:LGSSource}) = source.source

@inline _four_pupil_lgs_wavelength(source::LGSSource) = wavelength(source)
@inline _four_pupil_lgs_wavelength(
    source::FourPupilSpectralComponent{<:LGSSource}) = wavelength(source)

function prepare_four_pupil_lgs(source, input, front_end)
    return prepare_four_pupil_lgs(_four_pupil_lgs_source(source), source,
        input, front_end)
end

@inline prepare_four_pupil_lgs(::Nothing, source, input, front_end) =
    NoPreparedFourPupilLGS()

function prepare_four_pupil_lgs(source::LGSSource, spectral_source, input,
    front_end)
    return prepare_four_pupil_lgs(lgs_profile(source), source,
        _four_pupil_lgs_wavelength(spectral_source), input, front_end)
end

function prepare_four_pupil_lgs(::LGSProfileNone, source::LGSSource,
    wavelength_m, input, front_end)
    T = eltype(front_end.propagation.intensity)
    factor = T(lgs_elongation_factor(source))
    factor <= one(T) && return NoPreparedFourPupilLGS()
    sigma = T(0.5) * (factor - one(T))
    sigma <= zero(T) && return NoPreparedFourPupilLGS()
    half = max(1, ceil(Int, 2 * sigma))
    needed = 2 * half + 1
    host_kernel = Vector{T}(undef, needed)
    @inbounds for offset in -half:half
        host_kernel[offset + half + 1] =
            exp(-T(0.5) * (T(offset) / sigma)^2)
    end
    host_kernel ./= sum(host_kernel)
    kernel = similar(front_end.propagation.elongation_kernel, T, needed)
    copyto!(kernel, host_kernel)
    return PreparedFourPupilElongation(kernel, half)
end

function prepare_four_pupil_lgs(::LGSProfileNaProfile, source::LGSSource,
    wavelength_m, input, front_end)
    metadata = input.metadata
    resolution = metadata.dimensions[1]
    metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "LGS pupil input must be square"))
    metadata.sampling[1] == metadata.sampling[2] ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "LGS pupil input requires equal metric sampling on both axes"))
    propagation = front_end.propagation
    pad = size(propagation.field, 1)
    pupil_diameter = metadata.sampling[1] * resolution
    padding = pad / resolution
    pixel_scale = lgs_pixel_scale(pupil_diameter, padding, wavelength_m)
    kernel_fft = lgs_average_kernel_fft(pupil_diameter, source, pad,
        front_end.pupil_samples, pixel_scale, propagation.focal_field,
        propagation.fft_plan)
    return PreparedFourPupilSodiumProfile(kernel_fft)
end

@inline function apply_prepared_four_pupil_lgs!(
    ::NoPreparedFourPupilLGS, intensity, scratch, fft_buffer, fft_plan,
    ifft_buffer, ifft_plan)
    return intensity
end

function apply_prepared_four_pupil_lgs!(
    model::PreparedFourPupilElongation, intensity::AbstractMatrix{T},
    scratch::AbstractMatrix{T}, fft_buffer, fft_plan, ifft_buffer,
    ifft_plan) where {T<:AbstractFloat}
    n1, n2 = size(intensity)
    _apply_elongation!(execution_style(intensity), intensity, scratch,
        model.kernel, model.half_width, n1, n2)
    copyto!(intensity, scratch)
    return intensity
end

function apply_prepared_four_pupil_lgs!(
    model::PreparedFourPupilSodiumProfile, intensity::AbstractMatrix{T},
    scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan) where {
    T<:AbstractFloat,
}
    apply_lgs_convolution!(intensity, model.kernel_fft, fft_buffer, fft_plan,
        ifft_buffer, ifft_plan)
    return intensity
end

@inline four_pupil_bundle_input(input, ::Int) = input
@inline four_pupil_bundle_input(input::Union{Tuple,AbstractVector}, index::Int) =
    @inbounds input[index]

@inline form_four_pupil_bundle!(output, input, ::Tuple{}) = output

@inline function form_four_pupil_bundle!(output, input,
    plans::Tuple{P,Vararg{Any,N}}) where {P,N}
    index = length(output) - N
    component_input = four_pupil_bundle_input(input, index)
    form_wfs_optical_products!(output[index], component_input, first(plans))
    return form_four_pupil_bundle!(output, input, Base.tail(plans))
end

@inline four_pupil_path_sources(source::Asterism) = source.sources
@inline four_pupil_path_sources(source::ExtendedSource) =
    extended_source_asterism(source).sources

#
# WFS-specific focal-plane modulation policy and front-end integration
#
# Reusable modulation models and prepared optical quadrature are owned by
# `Optics`; this file retains WFS calibration and optics policy.
#

function legacy_modulation_policy(modulation::T,
    modulation_points::Union{Int,Nothing}, extra_modulation_factor::Int,
    modulation_phase_offset_rad::T,
    user_modulation_path) where {T<:AbstractFloat}
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
    return CircularModulation(modulation, samples,
        modulation_phase_offset_rad)
end

@inline calibration_modulation_policy(policy::SampledModulation,
    radius, phase_offset_rad) = policy

function calibration_modulation_policy(
    policy::Union{NoModulation,CircularModulation}, radius::T,
    phase_offset_rad::T) where {T<:AbstractFloat}
    iszero(radius) && return NoModulation()
    return CircularModulation(radius, modulation_point_count(policy),
        phase_offset_rad)
end

@inline _modulated_input_storage(input::PupilFunction) = input.opd
@inline _modulated_input_storage(input::ElectricField) = input.values

function require_modulated_wfs_input(input::PupilFunction)
    validate_wfs_optical_input(input)
    input.metadata.coordinate_domain isa MetricCoordinates ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "modulated WFS pupil input must use metric coordinates"))
    input.metadata.normalization isa DimensionlessNormalization ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "PupilFunction amplitude must be dimensionless"))
    input.metadata.spatial_measure isa PointSampledMeasure ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "PupilFunction amplitude must be point sampled"))
    input.metadata.coherence isa CoherentFieldCombination ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "modulated WFS pupil input must be coherent"))
    input.metadata.spectral isa AchromaticSpectralCoordinate ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "PupilFunction input must be achromatic"))
    return input
end

function require_modulated_wfs_input(input::ElectricField)
    validate_wfs_optical_input(input)
    input.metadata.coordinate_domain isa MetricCoordinates ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "modulated WFS electric-field input must use metric coordinates"))
    input.metadata.normalization isa PhotonRateNormalization ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "electric-field input must carry photon-rate normalization"))
    input.metadata.spatial_measure isa CellIntegratedMeasure ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "electric-field input must carry cell-integrated photon rate"))
    input.metadata.coherence isa CoherentFieldCombination ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "modulated WFS electric-field input must be coherent"))
    input.metadata.spectral isa MonochromaticChannel ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "electric-field input must declare one monochromatic channel"))
    return input
end

@inline function modulated_input_wavelength(input::ElectricField)
    return input.metadata.spectral.wavelength_m
end

function modulated_input_wavelength(input::PupilFunction, source)
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "WFS optics require an illumination source for PupilFunction input"))
    return wavelength(source)
end

function four_pupil_propagation_workspace end

@inline modulated_wfs_propagation_storage(front_end) =
    front_end.propagation.field

function require_modulated_wfs_domains(front_end, input, output::IntensityMap)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:wfs_optics, :backend,
            "modulated WFS input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:wfs_optics, :device,
            "modulated WFS input and output occupy different devices"))
    storage = modulated_wfs_propagation_storage(front_end)
    typeof(input.metadata.backend) === typeof(backend(storage)) ||
        throw(WFSPreparationError(:wfs_optics, :backend,
            "modulated WFS input and propagation backends differ"))
    input.metadata.device == compute_device(storage) ||
        throw(WFSPreparationError(:wfs_optics, :device,
            "modulated WFS input and propagation occupy different devices"))
    return nothing
end

function require_four_pupil_rate_map(output::IntensityMap, expected_size,
    wavelength_m)
    validate_wfs_optical_products(output)
    output.metadata.coordinate_domain isa NormalizedPupilCoordinates ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "four-pupil detector output must use normalized pupil coordinates"))
    output.metadata.spatial_measure isa CellIntegratedMeasure ||
        throw(WFSPreparationError(:wfs_optics, :radiometry,
            "four-pupil detector output must carry cell-integrated rate"))
    size(output.values) == expected_size || throw(WFSPreparationError(
        :wfs_optics, :shape,
        "four-pupil detector output has the wrong prepared dimensions"))
    channel = output.metadata.spectral
    channel isa MonochromaticChannel &&
        channel.wavelength_m == wavelength_m ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
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

"""
    AbstractPreparedFourPupilLGS

Internal prepared interface for four-pupil LGS image formation. Implementations
define `apply_prepared_four_pupil_lgs!` and
`_require_exact_prepared_four_pupil_lgs_target`. They own only run-immutable
execution data. Intensity products and FFT scratch remain distinct,
caller-owned, single-writer storage on the exact prepared target. A prepared
value is reentrant only across distinct workspaces. Preparation and target
validation raise `WFSPreparationError` before optical-product mutation.
"""
abstract type AbstractPreparedFourPupilLGS end

struct NoPreparedFourPupilLGS <: AbstractPreparedFourPupilLGS end

struct PreparedFourPupilElongation{K<:AbstractVector} <:
    AbstractPreparedFourPupilLGS
    kernel::K
    half_width::Int
end

struct PreparedFourPupilSodiumLayerProfile{K<:AbstractMatrix} <:
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
    return prepare_four_pupil_lgs(sodium_layer_profile_style(source), source,
        _four_pupil_lgs_wavelength(spectral_source), input, front_end)
end

function prepare_four_pupil_lgs(::NoSodiumLayerProfileStyle, source::LGSSource,
    wavelength_m, input, front_end)
    propagation = four_pupil_propagation_workspace(front_end)
    T = eltype(propagation.intensity)
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
    kernel = similar(propagation.elongation_kernel, T, needed)
    copyto!(kernel, host_kernel)
    return PreparedFourPupilElongation(kernel, half)
end

function prepare_four_pupil_lgs(::SampledSodiumLayerProfileStyle, source::LGSSource,
    wavelength_m, input, front_end)
    metadata = input.metadata
    resolution = metadata.dimensions[1]
    metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "LGS pupil input must be square"))
    metadata.sampling[1] == metadata.sampling[2] ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "LGS pupil input requires equal metric sampling on both axes"))
    propagation = four_pupil_propagation_workspace(front_end)
    pad = size(propagation.field, 1)
    pupil_diameter = metadata.sampling[1] * resolution
    padding = pad / resolution
    pixel_scale = lgs_pixel_scale(pupil_diameter, padding, wavelength_m)
    kernel_fft = lgs_average_kernel_fft(pupil_diameter, source, pad,
        front_end.pupil_samples, pixel_scale, propagation.focal_field,
        propagation.fft_plan)
    return PreparedFourPupilSodiumLayerProfile(kernel_fft)
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
    model::PreparedFourPupilSodiumLayerProfile, intensity::AbstractMatrix{T},
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

@inline validate_four_pupil_bundle_binding(output, input, ::Tuple{}) = nothing

@inline function validate_four_pupil_bundle_binding(output, input,
    plans::Tuple{P,Vararg{Any,N}}) where {P,N}
    index = length(output) - N
    component_input = four_pupil_bundle_input(input, index)
    validate_wfs_optics_binding(output[index], component_input,
        first(plans))
    return validate_four_pupil_bundle_binding(output, input,
        Base.tail(plans))
end

@inline four_pupil_path_sources(source::Asterism) = source.sources
@inline four_pupil_path_sources(source::ExtendedSource) =
    extended_source_asterism(source).sources

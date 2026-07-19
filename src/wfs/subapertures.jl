abstract type AbstractSlopeExtractionModel end
abstract type AbstractValidSubaperturePolicy end

struct GeometryValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    threshold::T
end

struct FluxThresholdValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    light_ratio::T
end

function GeometryValidSubapertures(; threshold::Real=0.1,
    T::Type{<:AbstractFloat}=typeof(float(threshold)))
    value = T(threshold)
    isfinite(value) && zero(T) <= value <= one(T) ||
        throw(InvalidConfiguration(
            "GeometryValidSubapertures threshold must lie in [0, 1]"))
    return GeometryValidSubapertures{T}(value)
end

function FluxThresholdValidSubapertures(; light_ratio::Real=0.5,
    T::Type{<:AbstractFloat}=typeof(float(light_ratio)))
    value = T(light_ratio)
    isfinite(value) && zero(T) <= value <= one(T) ||
        throw(InvalidConfiguration(
            "FluxThresholdValidSubapertures light_ratio must lie in [0, 1]"))
    return FluxThresholdValidSubapertures{T}(value)
end

struct CenterOfGravityExtraction{T<:AbstractFloat,W<:Union{Nothing,AbstractMatrix{T}}} <: AbstractSlopeExtractionModel
    threshold::T
    window::W
end

mutable struct SubapertureLayoutState
    revision::UInt
end

function CenterOfGravityExtraction(threshold::Real; window=nothing,
    T::Type{<:AbstractFloat}=typeof(float(threshold)))
    value = T(threshold)
    isfinite(value) && zero(T) <= value <= one(T) ||
        throw(InvalidConfiguration(
            "CenterOfGravityExtraction threshold must lie in [0, 1]"))
    window === nothing || throw(UnsupportedAlgorithm(
        "windowed center-of-gravity extraction is not implemented"))
    return CenterOfGravityExtraction{T,Nothing}(value, nothing)
end

"""
Fixed subaperture geometry with caller-owned mutable mask storage.

The execution mask is authoritative at construction. Its host mirror and
column-major `CartesianIndex` list are synchronized before the layout is
returned. Maintained updates synchronize all three representations and advance
a cold-configuration revision; direct storage mutation is unsupported.
"""
struct SubapertureLayout{T<:AbstractFloat,A<:AbstractMatrix{Bool},M<:Matrix{Bool},V<:Vector{CartesianIndex{2}}}
    n_subap::Int
    subap_pixels::Int
    pitch_m::T
    threshold::T
    valid_mask::A
    valid_mask_host::M
    valid_indices_host::V
    state::SubapertureLayoutState
end

mutable struct SubapertureCalibration{T<:AbstractFloat,
    R<:AbstractMatrix{T},V<:AbstractVector{T},
    E<:AbstractSlopeExtractionModel,U}
    extraction::E
    reference_signal_2d::R
    reference_signal_host::V
    centroid_response::T
    output_units::U
    calibrated::Bool
    wavelength::T
    signature::UInt
    revision::UInt
end

function SubapertureLayout(n_subap::Int, pupil_resolution::Int, diameter::Real, threshold::Real,
    valid_mask::AbstractMatrix{Bool}, valid_mask_host::Matrix{Bool})
    n_subap > 0 || throw(InvalidConfiguration(
        "SubapertureLayout n_subap must be positive"))
    pupil_resolution > 0 || throw(InvalidConfiguration(
        "SubapertureLayout pupil_resolution must be positive"))
    pupil_resolution % n_subap == 0 || throw(InvalidConfiguration(
        "SubapertureLayout pupil_resolution must be divisible by n_subap"))
    size(valid_mask) == (n_subap, n_subap) || throw(DimensionMismatchError(
        "SubapertureLayout valid_mask must have size (n_subap, n_subap)"))
    size(valid_mask_host) == size(valid_mask) || throw(DimensionMismatchError(
        "SubapertureLayout host and execution masks must have the same size"))
    T = promote_type(typeof(float(diameter)), typeof(float(threshold)))
    diameter_m = T(diameter)
    threshold_fraction = T(threshold)
    isfinite(diameter_m) && diameter_m > zero(T) ||
        throw(InvalidConfiguration(
            "SubapertureLayout diameter must be finite and positive"))
    isfinite(threshold_fraction) &&
        zero(T) <= threshold_fraction <= one(T) ||
        throw(InvalidConfiguration(
            "SubapertureLayout threshold must lie in [0, 1]"))
    subap_pixels = div(pupil_resolution, n_subap)
    pitch_m = diameter_m / T(n_subap)
    valid_indices_host = CartesianIndex{2}[]
    sizehint!(valid_indices_host, length(valid_mask_host))
    layout = SubapertureLayout{T, typeof(valid_mask), typeof(valid_mask_host), Vector{CartesianIndex{2}}}(
        n_subap,
        subap_pixels,
        pitch_m,
        threshold_fraction,
        valid_mask,
        valid_mask_host,
        valid_indices_host,
        SubapertureLayoutState(UInt(0)),
    )
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    return layout
end

@inline subaperture_layout_revision(layout::SubapertureLayout) =
    layout.state.revision

@inline function _advance_subaperture_layout_revision!(
    layout::SubapertureLayout)
    layout.state.revision += UInt(1)
    return layout.state.revision
end

@inline function _copy_valid_mask_to_host!(host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    return _copy_valid_mask_to_host!(execution_style(mask), host, mask)
end

@inline function _copy_valid_mask_to_host!(::ScalarCPUStyle, host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    copyto!(host, mask)
    return host
end

@inline function _copy_valid_mask_to_host!(::ExecutionStyle, host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    copyto!(host, Array(mask))
    return host
end

@inline function _refresh_valid_indices_host!(layout::SubapertureLayout)
    indices = layout.valid_indices_host
    resize!(indices, 0)
    @inbounds for I in CartesianIndices(layout.valid_mask_host)
        layout.valid_mask_host[I] && push!(indices, I)
    end
    return indices
end

function SubapertureCalibration(reference_signal_2d::AbstractMatrix{T},
    reference_signal_host::AbstractVector{T},
    extraction::AbstractSlopeExtractionModel=CenterOfGravityExtraction(
        T(0.01); T=T);
    output_units=:pixel) where {T<:AbstractFloat}
    size(reference_signal_2d, 2) == 2 ||
        throw(DimensionMismatchError(
            "SubapertureCalibration reference must have two component columns"))
    n_lenslet_values = size(reference_signal_2d, 1)
    n_subap = isqrt(n_lenslet_values)
    n_subap * n_subap == n_lenslet_values ||
        throw(DimensionMismatchError(
            "SubapertureCalibration reference rows must form a square lenslet grid"))
    length(reference_signal_host) == length(reference_signal_2d) ||
        throw(DimensionMismatchError(
            "SubapertureCalibration host reference length must match reference storage"))
    _require_declared_wfs_units(output_units, :estimation)
    reference_host = vec(Array(reference_signal_2d))
    all(isfinite, reference_host) || throw(InvalidConfiguration(
        "SubapertureCalibration reference must contain only finite values"))
    copyto!(reference_signal_host, reference_host)
    return SubapertureCalibration{T,typeof(reference_signal_2d),
        typeof(reference_signal_host),typeof(extraction),typeof(output_units)}(
        extraction,
        reference_signal_2d,
        reference_signal_host,
        one(T),
        output_units,
        false,
        zero(T),
        UInt(0),
        UInt(0),
    )
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool})
    build_mask!(layout.valid_mask, SubapertureGridMask(threshold=layout.threshold, T=typeof(layout.threshold)), pupil)
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool},
    policy::GeometryValidSubapertures)
    threshold = convert(typeof(layout.threshold), policy.threshold)
    build_mask!(layout.valid_mask,
        SubapertureGridMask(threshold=threshold,
            T=typeof(layout.threshold)), pupil)
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout,
    support_map::AbstractMatrix{T},
    policy::GeometryValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(support_map) == (n_sub * sub, n_sub * sub) ||
        throw(DimensionMismatchError(
            "support map size must match subaperture layout"))
    support_host = _host_support_map(execution_style(support_map),
        support_map)
    threshold = convert(Float64, policy.threshold)
    denominator = sub * sub
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        illuminated = count(!iszero,
            @view support_host[xs:xe, ys:ye])
        layout.valid_mask_host[i, j] =
            illuminated / denominator >= threshold
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

@inline _host_support_map(::ScalarCPUStyle, support_map::AbstractMatrix) = support_map
@inline _host_support_map(::ExecutionStyle, support_map::AbstractMatrix) = Array(support_map)

function update_subaperture_layout!(layout::SubapertureLayout, support_map::AbstractMatrix{T},
    policy::FluxThresholdValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(support_map, 1) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    size(support_map, 2) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    support_host = _host_support_map(execution_style(support_map), support_map)
    peak = zero(eltype(support_host))
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        total = sum(@view support_host[xs:xe, ys:ye])
        peak = max(peak, total)
    end
    cutoff = convert(eltype(support_host), policy.light_ratio) * peak
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        layout.valid_mask_host[i, j] = sum(@view support_host[xs:xe, ys:ye]) >= cutoff
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout_from_amplitude!(
    layout::SubapertureLayout, amplitude::AbstractMatrix{T},
    policy::FluxThresholdValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(amplitude) == (n_sub * sub, n_sub * sub) ||
        throw(DimensionMismatchError(
            "pupil amplitude size must match subaperture layout"))
    amplitude_host = _host_support_map(execution_style(amplitude),
        amplitude)
    peak = zero(eltype(amplitude_host))
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        total = sum(abs2, @view amplitude_host[xs:xe, ys:ye])
        peak = max(peak, total)
    end
    cutoff = convert(eltype(amplitude_host), policy.light_ratio) * peak
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        layout.valid_mask_host[i, j] =
            sum(abs2, @view amplitude_host[xs:xe, ys:ye]) >= cutoff
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function set_reference_signal!(calibration::SubapertureCalibration, signal_2d::AbstractMatrix)
    size(signal_2d) == size(calibration.reference_signal_2d) ||
        throw(DimensionMismatchError("reference signal shape must match SubapertureCalibration storage"))
    copyto!(calibration.reference_signal_2d, signal_2d)
    copyto!(calibration.reference_signal_host, vec(Array(calibration.reference_signal_2d)))
    calibration.calibrated = false
    calibration.revision += UInt(1)
    return calibration
end

"""
    set_subaperture_calibration!(calibration, reference_signal;
        centroid_response, output_units=:pixel, wavelength, signature=0)

Validate and install an explicit Shack-Hartmann reference and centroid
calibration under the package's single-writer state contract.
`centroid_response` is the raw reference-subtracted centroid response per
reported output unit. Prepared estimators bind the resulting revision, output
units, and calibration state and reject subsequent changes made through this
API.
"""
function set_subaperture_calibration!(
    calibration::SubapertureCalibration,
    reference_signal::AbstractMatrix;
    centroid_response::Real,
    output_units=calibration.output_units,
    wavelength::Real,
    signature::Integer=0)
    size(reference_signal) == size(calibration.reference_signal_2d) ||
        throw(DimensionMismatchError(
            "reference signal shape must match SubapertureCalibration storage"))
    T = typeof(calibration.centroid_response)
    response = T(centroid_response)
    isfinite(response) && response != zero(T) ||
        throw(InvalidConfiguration(
            "centroid_response must be finite and nonzero"))
    wavelength_m = T(wavelength)
    isfinite(wavelength_m) && wavelength_m > zero(T) ||
        throw(InvalidConfiguration(
            "calibration wavelength must be finite and positive"))
    reference_host = Array(reference_signal)
    all(isfinite, reference_host) || throw(InvalidConfiguration(
        "reference signal must contain only finite values"))
    _require_declared_wfs_units(output_units, :estimation)
    typeof(output_units) === typeof(calibration.output_units) ||
        throw(InvalidConfiguration(
            "calibration output-unit descriptor type cannot change in place"))
    signature >= 0 || throw(InvalidConfiguration(
        "calibration signature must be nonnegative"))
    signature <= typemax(UInt) || throw(InvalidConfiguration(
        "calibration signature exceeds UInt range"))
    signature_value = UInt(signature)

    copyto!(calibration.reference_signal_2d, reference_host)
    copyto!(calibration.reference_signal_host, vec(reference_host))
    calibration.centroid_response = response
    calibration.output_units = output_units
    calibration.calibrated = true
    calibration.wavelength = wavelength_m
    calibration.signature = signature_value
    calibration.revision += UInt(1)
    return calibration
end

@inline centroid_cutoff(model::CenterOfGravityExtraction{T}, peak::T) where {T<:AbstractFloat} =
    peak <= zero(T) ? zero(T) : model.threshold * peak

@inline slope_extraction_model(calibration::SubapertureCalibration) = calibration.extraction
@inline reference_signal(calibration::SubapertureCalibration) = calibration.reference_signal_2d
@inline reference_signal_host(calibration::SubapertureCalibration) = calibration.reference_signal_host
@inline valid_subaperture_indices(layout::SubapertureLayout) = layout.valid_indices_host
@inline n_valid_subapertures(layout::SubapertureLayout) = length(layout.valid_indices_host)

function subaperture_layout end
function subaperture_calibration end

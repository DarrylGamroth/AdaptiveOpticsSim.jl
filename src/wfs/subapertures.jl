abstract type AbstractSlopeExtractionModel end
abstract type AbstractValidSubaperturePolicy end

struct GeometryValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    threshold::T
end

struct FluxThresholdValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    light_ratio::T
end

GeometryValidSubapertures(; threshold::Real=0.1, T::Type{<:AbstractFloat}=typeof(float(threshold))) =
    GeometryValidSubapertures{T}(T(threshold))

FluxThresholdValidSubapertures(; light_ratio::Real=0.5, T::Type{<:AbstractFloat}=typeof(float(light_ratio))) =
    FluxThresholdValidSubapertures{T}(T(light_ratio))

struct CenterOfGravityExtraction{T<:AbstractFloat,W<:Union{Nothing,AbstractMatrix{T}}} <: AbstractSlopeExtractionModel
    threshold::T
    window::W
end

CenterOfGravityExtraction(threshold::Real; window=nothing, T::Type{<:AbstractFloat}=typeof(float(threshold))) =
    CenterOfGravityExtraction{T, typeof(window)}(T(threshold), window)

mutable struct SubapertureLayout{T<:AbstractFloat,A<:AbstractMatrix{Bool},M<:Matrix{Bool},V<:Vector{CartesianIndex{2}}}
    n_subap::Int
    subap_pixels::Int
    pitch_m::T
    threshold::T
    valid_mask::A
    valid_mask_host::M
    valid_indices_host::V
end

mutable struct SubapertureCalibration{T<:AbstractFloat,R<:AbstractMatrix{T},V<:AbstractVector{T},E<:AbstractSlopeExtractionModel}
    extraction::E
    reference_signal_2d::R
    reference_signal_host::V
    slopes_units::T
    calibrated::Bool
    wavelength::T
    signature::UInt
end

function SubapertureLayout(n_subap::Int, pupil_resolution::Int, diameter::Real, threshold::Real,
    valid_mask::AbstractMatrix{Bool}, valid_mask_host::Matrix{Bool})
    T = promote_type(typeof(diameter), typeof(threshold))
    subap_pixels = div(pupil_resolution, n_subap)
    pitch_m = T(diameter / n_subap)
    return SubapertureLayout{T, typeof(valid_mask), typeof(valid_mask_host), Vector{CartesianIndex{2}}}(
        n_subap,
        subap_pixels,
        pitch_m,
        T(threshold),
        valid_mask,
        valid_mask_host,
        CartesianIndex{2}[],
    )
end

function SubapertureCalibration(reference_signal_2d::AbstractMatrix{T}, reference_signal_host::AbstractVector{T},
    extraction::AbstractSlopeExtractionModel=CenterOfGravityExtraction(T(0.01); T=T)) where {T<:AbstractFloat}
    return SubapertureCalibration{T, typeof(reference_signal_2d), typeof(reference_signal_host), typeof(extraction)}(
        extraction,
        reference_signal_2d,
        reference_signal_host,
        one(T),
        false,
        zero(T),
        UInt(0),
    )
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool})
    build_mask!(layout.valid_mask, SubapertureGridMask(threshold=layout.threshold, T=typeof(layout.threshold)), pupil)
    copyto!(layout.valid_mask_host, Array(layout.valid_mask))
    resize!(layout.valid_indices_host, 0)
    append!(layout.valid_indices_host, findall(layout.valid_mask_host))
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool},
    policy::GeometryValidSubapertures)
    layout.threshold = convert(typeof(layout.threshold), policy.threshold)
    build_mask!(layout.valid_mask, SubapertureGridMask(threshold=layout.threshold, T=typeof(layout.threshold)), pupil)
    copyto!(layout.valid_mask_host, Array(layout.valid_mask))
    resize!(layout.valid_indices_host, 0)
    append!(layout.valid_indices_host, findall(layout.valid_mask_host))
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout, support_map::AbstractMatrix{T},
    policy::FluxThresholdValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(support_map, 1) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    size(support_map, 2) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    support_host = Array(support_map)
    flux_per_subap = Matrix{Float64}(undef, n_sub, n_sub)
    peak = 0.0
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        total = sum(@view support_host[xs:xe, ys:ye])
        flux_per_subap[i, j] = total
        peak = max(peak, total)
    end
    cutoff = convert(Float64, policy.light_ratio) * peak
    @inbounds for i in 1:n_sub, j in 1:n_sub
        layout.valid_mask_host[i, j] = flux_per_subap[i, j] >= cutoff
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    resize!(layout.valid_indices_host, 0)
    append!(layout.valid_indices_host, findall(layout.valid_mask_host))
    return layout
end

function set_reference_signal!(calibration::SubapertureCalibration, signal_2d::AbstractMatrix)
    size(signal_2d) == size(calibration.reference_signal_2d) ||
        throw(DimensionMismatchError("reference signal shape must match SubapertureCalibration storage"))
    copyto!(calibration.reference_signal_2d, signal_2d)
    copyto!(calibration.reference_signal_host, vec(Array(calibration.reference_signal_2d)))
    return calibration
end

function set_calibration_state!(calibration::SubapertureCalibration;
    slopes_units::Real=calibration.slopes_units,
    calibrated::Bool=calibration.calibrated,
    wavelength::Real=calibration.wavelength,
    signature::UInt=calibration.signature)
    T = typeof(calibration.slopes_units)
    calibration.slopes_units = T(slopes_units)
    calibration.calibrated = calibrated
    calibration.wavelength = T(wavelength)
    calibration.signature = signature
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

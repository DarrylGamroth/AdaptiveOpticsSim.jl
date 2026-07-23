#
# Shack-Hartmann wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("shack_hartmann/setup.jl")
include("shack_hartmann/measure.jl")
include("shack_hartmann/stacks.jl")
include("shack_hartmann/signals.jl")
include("shack_hartmann/lgs.jl")
include("shack_hartmann/stages.jl")

@inline slopes(wfs::ShackHartmannWFS) = wfs.estimator.slopes
@inline valid_subaperture_mask(wfs::ShackHartmannWFS) =
    wfs.front_end.layout.valid_mask
@inline reference_signal(wfs::ShackHartmannWFS) = wfs.calibration.reference_signal_2d

@kernel function shack_hartmann_detector_image_kernel!(image, spot_cube,
    n_sub::Int, n_axis_1::Int, n_axis_2::Int, gap::Int, gap_value)
    axis_1_index, axis_2_index = @index(Global, NTuple)
    pitch_1 = n_axis_1 + gap
    pitch_2 = n_axis_2 + gap
    lenslet_i = (axis_1_index - 1) ÷ pitch_1 + 1
    lenslet_j = (axis_2_index - 1) ÷ pitch_2 + 1
    local_i = (axis_1_index - 1) % pitch_1 + 1
    local_j = (axis_2_index - 1) % pitch_2 + 1
    if lenslet_i <= n_sub && lenslet_j <= n_sub &&
            local_i <= n_axis_1 && local_j <= n_axis_2
        idx = sh_lenslet_index(lenslet_i, lenslet_j, n_sub)
        @inbounds image[axis_1_index, axis_2_index] =
            detector_output_value(eltype(image),
                spot_cube[idx, local_i, local_j])
    else
        @inbounds image[axis_1_index, axis_2_index] =
            detector_output_value(eltype(image), gap_value)
    end
end

@kernel function shack_hartmann_detector_image_copy_kernel!(image, spot_cube,
    n_sub::Int, n_axis_1::Int, n_axis_2::Int, gap::Int)
    lenslet_i, lenslet_j, local_i, local_j = @index(Global, NTuple)
    idx = sh_lenslet_index(lenslet_i, lenslet_j, n_sub)
    axis_1_index = (lenslet_i - 1) * (n_axis_1 + gap) + local_i
    axis_2_index = (lenslet_j - 1) * (n_axis_2 + gap) + local_j
    @inbounds image[axis_1_index, axis_2_index] =
        spot_cube[idx, local_i, local_j]
end

"""
    shack_hartmann_detector_image(spot_cube, n_lenslets; gap=0, gap_value=0)
    shack_hartmann_detector_image(wfs; gap=0, gap_value=0)

Tile a Shack-Hartmann lenslet spot cube into a 2-D detector image.
"""
function shack_hartmann_detector_image(spot_cube::AbstractArray, n_lenslets::Integer; gap::Integer=0, gap_value=0,
    output_type::Union{Nothing,DataType}=nothing)
    n_sub = Int(n_lenslets)
    n_sub > 0 || throw(ArgumentError("n_lenslets must be positive"))
    ndims(spot_cube) == 3 ||
        throw(DimensionMismatch("spot cube must have shape (n_lenslets^2, n_axis_1, n_axis_2)"))
    n_spots, n_axis_1, n_axis_2 = size(spot_cube)
    n_spots == n_sub * n_sub ||
        throw(DimensionMismatch("spot cube first dimension must equal n_lenslets^2; got $n_spots for n_lenslets=$n_sub"))
    gap_px = Int(gap)
    gap_px >= 0 || throw(ArgumentError("gap must be non-negative"))
    T = output_type === nothing ? promote_type(eltype(spot_cube), typeof(gap_value)) : output_type
    image = similar(spot_cube, T,
        n_sub * n_axis_1 + (n_sub - 1) * gap_px,
        n_sub * n_axis_2 + (n_sub - 1) * gap_px)
    return shack_hartmann_detector_image!(image, spot_cube, n_sub; gap=gap_px, gap_value=gap_value)
end

function shack_hartmann_detector_image!(image::AbstractMatrix, spot_cube::AbstractArray, n_lenslets::Integer;
    gap::Integer=0, gap_value=0)
    image_style = execution_style(image)
    spot_style = execution_style(spot_cube)
    n_sub = Int(n_lenslets)
    n_sub > 0 || throw(ArgumentError("n_lenslets must be positive"))
    ndims(spot_cube) == 3 ||
        throw(DimensionMismatch("spot cube must have shape (n_lenslets^2, n_axis_1, n_axis_2)"))
    n_spots, n_axis_1, n_axis_2 = size(spot_cube)
    n_spots == n_sub * n_sub ||
        throw(DimensionMismatch("spot cube first dimension must equal n_lenslets^2; got $n_spots for n_lenslets=$n_sub"))
    gap_px = Int(gap)
    gap_px >= 0 || throw(ArgumentError("gap must be non-negative"))
    size(image) == (
        n_sub * n_axis_1 + (n_sub - 1) * gap_px,
        n_sub * n_axis_2 + (n_sub - 1) * gap_px,
    ) ||
        throw(DimensionMismatch("output image has size $(size(image)); expected " *
            "$((n_sub * n_axis_1 + (n_sub - 1) * gap_px, n_sub * n_axis_2 + (n_sub - 1) * gap_px))"))

    _shack_hartmann_detector_image!(image_style, spot_style, image,
        spot_cube, n_sub, n_axis_1, n_axis_2, gap_px, gap_value)
    return image
end

function _shack_hartmann_detector_image!(::ExecutionStyle, ::ExecutionStyle, image::AbstractMatrix,
    spot_cube::AbstractArray, n_sub::Int, n_axis_1::Int, n_axis_2::Int,
    gap::Int, gap_value)
    throw(InvalidConfiguration("Shack-Hartmann detector image output and spot cube must use the same execution backend"))
end

function _shack_hartmann_detector_image!(::ScalarCPUStyle, ::ScalarCPUStyle, image::AbstractMatrix,
    spot_cube::AbstractArray, n_sub::Int, n_axis_1::Int, n_axis_2::Int,
    gap::Int, gap_value)
    T = eltype(image)
    fill!(image, detector_output_value(T, gap_value))
    @inbounds for j in 1:n_sub, i in 1:n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        axis_1_start = (i - 1) * (n_axis_1 + gap) + 1
        axis_2_start = (j - 1) * (n_axis_2 + gap) + 1
        for local_j in 1:n_axis_2, local_i in 1:n_axis_1
            image[axis_1_start + local_i - 1,
                axis_2_start + local_j - 1] = detector_output_value(T,
                    spot_cube[idx, local_i, local_j])
        end
    end
    return image
end

function _shack_hartmann_detector_image!(style::AcceleratorStyle{B}, ::AcceleratorStyle{B},
    image::AbstractMatrix, spot_cube::AbstractArray, n_sub::Int,
    n_axis_1::Int, n_axis_2::Int, gap::Int, gap_value) where {B}
    output_gap_value = eltype(image)(gap_value)
    if eltype(image) === eltype(spot_cube)
        gap > 0 && fill!(image, output_gap_value)
        launch_kernel!(style, shack_hartmann_detector_image_copy_kernel!, image, spot_cube,
            n_sub, n_axis_1, n_axis_2, gap;
            ndrange=(n_sub, n_sub, n_axis_1, n_axis_2))
    else
        launch_kernel!(style, shack_hartmann_detector_image_kernel!, image, spot_cube,
            n_sub, n_axis_1, n_axis_2, gap, output_gap_value;
            ndrange=size(image))
    end
    return image
end

@inline shack_hartmann_detector_image(wfs::ShackHartmannWFS; kwargs...) =
    shack_hartmann_detector_image(sh_exported_spot_cube(wfs), n_lenslets(wfs); kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, ::Nothing; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, det::AbstractDetector; kwargs...) =
    shack_hartmann_detector_image(wfs; output_type=detector_output_type(det), kwargs...)

@inline wfs_output_frame(wfs::ShackHartmannWFS{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame(wfs::ShackHartmannWFS{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmannWFS{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmannWFS{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline function wfs_output_metadata(wfs::ShackHartmannWFS)
    layout = wfs.front_end.layout
    calibration = subaperture_calibration(wfs)
    return (
        n_lenslets=layout.n_subap,
        n_valid_subap=n_valid_subapertures(layout),
        subap_pixels=layout.subap_pixels,
        pitch_m=layout.pitch_m,
        centroid_response=calibration.centroid_response,
        calibrated=calibration.calibrated,
    )
end

@inline supports_prepared_runtime(::ShackHartmannWFS{<:Diffractive}, ::AbstractSource) = true
@inline supports_prepared_runtime(wfs::ShackHartmannWFS{<:Diffractive},
    src::SpectralSource) = sh_has_common_spectral_grid(wfs, src)
@inline supports_prepared_runtime(::ShackHartmannWFS{<:Diffractive}, ::Asterism) = true
@inline supports_detector_output(::ShackHartmannWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::ShackHartmannWFS, ::Asterism) = true
@inline supports_stacked_sources(::ShackHartmannWFS, ::SpectralSource) = true
@inline supports_stacked_sources(wfs::ShackHartmannWFS{<:Diffractive},
    src::SpectralSource) = sh_has_common_spectral_grid(wfs, src)
@inline supports_stacked_sources(::ShackHartmannWFS, ::ExtendedSource) = true
@inline supports_grouped_execution(::ShackHartmannWFS{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(wfs::ShackHartmannWFS{<:Diffractive},
    src::SpectralSource) = sh_has_common_spectral_grid(wfs, src)
@inline supports_grouped_execution(::ShackHartmannWFS{<:Diffractive}, ::ExtendedSource) = true

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, pupil::PupilFunction, src::AbstractSource)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, pupil::PupilFunction, src::SpectralSource)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, pupil::PupilFunction, ast::Asterism)
    common_source = common_wfs_calibration_source(ast, "ShackHartmannWFS")
    prepare_sampling!(wfs, pupil, common_source)
    ensure_sh_calibration!(wfs, pupil, common_source)
    return wfs
end

@inline function _measure_for_calibration!(wfs::ShackHartmannWFS{<:Diffractive}, pupil::PupilFunction, src::AbstractSource)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    return measure!(wfs, pupil, src)
end

@inline function _measure_for_calibration!(wfs::ShackHartmannWFS{<:Diffractive}, pupil::PupilFunction, src::SpectralSource)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    return measure!(wfs, pupil, src)
end

#
# Wavefront-sensor interface accessors
#
# This file defines the maintained subsystem-level accessors for `AbstractWFS`.
# The goal is to keep runtime, calibration, and tests off direct `state` field
# access when the accessed quantity is part of the supported WFS contract.
#

"""
    slopes(wfs::AbstractWFS)

Return the current exported 1-D signal vector for a wavefront sensor.

For geometric and centroid-based sensors this is the slope vector. For sensors
such as `ZernikeWFS` and `CurvatureWFS`, the maintained runtime contract still
uses the `slopes` name even though the stored values represent normalized
signal samples rather than geometric slopes.
"""
@inline slopes(wfs::AbstractWFS) = wfs.state.slopes

"""
    valid_subaperture_mask(wfs::AbstractWFS)

Return the maintained valid-subaperture/support mask when the WFS exposes one.
Returns `nothing` for sensor families that do not export a stable valid-mask
surface.
"""
@inline valid_subaperture_mask(::AbstractWFS) = nothing

"""
    reference_signal(wfs::AbstractWFS)

Return the stored reference signal surface used by the WFS calibration path, or
`nothing` when the family does not expose one as part of the maintained
interface.
"""
@inline reference_signal(::AbstractWFS) = nothing

"""
    camera_frame(wfs::AbstractWFS)

Return the maintained internal camera/readout frame for detector-like WFS
families, or `nothing` when the family does not expose a stable camera-frame
surface.
"""
@inline camera_frame(::AbstractWFS) = nothing

"""
    wfs_detector_image(wfs)
    wfs_detector_image(wfs, det)

Return the detector-image product for a wavefront sensor.

For WFS families with a native 2-D camera frame, this returns `camera_frame(wfs)`
when no detector is provided and `output_frame(det)` after detector-coupled
measurement when a detector is provided. Shack-Hartmann WFS objects store
detector-coupled lenslet spots as a cube, so this accessor returns a tiled 2-D
detector mosaic.
"""
function wfs_detector_image(wfs::AbstractWFS)
    frame = camera_frame(wfs)
    isnothing(frame) && throw(InvalidConfiguration("WFS of type $(typeof(wfs)) does not expose a detector image"))
    return frame
end

@inline wfs_detector_image(::AbstractWFS, det::AbstractDetector) = output_frame(det)

@kernel function shack_hartmann_detector_image_kernel!(image, spot_cube, n_sub::Int, n_y::Int, n_x::Int, gap::Int, gap_value)
    y, x = @index(Global, NTuple)
    pitch_y = n_y + gap
    pitch_x = n_x + gap
    tile_y = (y - 1) ÷ pitch_y + 1
    tile_x = (x - 1) ÷ pitch_x + 1
    local_y = (y - 1) % pitch_y + 1
    local_x = (x - 1) % pitch_x + 1
    if tile_y <= n_sub && tile_x <= n_sub && local_y <= n_y && local_x <= n_x
        idx = (tile_y - 1) * n_sub + tile_x
        @inbounds image[y, x] = detector_output_value(eltype(image), spot_cube[idx, local_y, local_x])
    else
        @inbounds image[y, x] = detector_output_value(eltype(image), gap_value)
    end
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
        throw(DimensionMismatch("spot cube must have shape (n_lenslets^2, n_pix_subap_y, n_pix_subap_x)"))
    n_spots, n_y, n_x = size(spot_cube)
    n_spots == n_sub * n_sub ||
        throw(DimensionMismatch("spot cube first dimension must equal n_lenslets^2; got $n_spots for n_lenslets=$n_sub"))
    gap_px = Int(gap)
    gap_px >= 0 || throw(ArgumentError("gap must be non-negative"))
    T = output_type === nothing ? promote_type(eltype(spot_cube), typeof(gap_value)) : output_type
    image = similar(spot_cube, T, n_sub * n_y + (n_sub - 1) * gap_px,
        n_sub * n_x + (n_sub - 1) * gap_px)
    return shack_hartmann_detector_image!(image, spot_cube, n_sub; gap=gap_px, gap_value=gap_value)
end

function shack_hartmann_detector_image!(image::AbstractMatrix, spot_cube::AbstractArray, n_lenslets::Integer;
    gap::Integer=0, gap_value=0)
    image_style = execution_style(image)
    spot_style = execution_style(spot_cube)
    n_sub = Int(n_lenslets)
    n_sub > 0 || throw(ArgumentError("n_lenslets must be positive"))
    ndims(spot_cube) == 3 ||
        throw(DimensionMismatch("spot cube must have shape (n_lenslets^2, n_pix_subap_y, n_pix_subap_x)"))
    n_spots, n_y, n_x = size(spot_cube)
    n_spots == n_sub * n_sub ||
        throw(DimensionMismatch("spot cube first dimension must equal n_lenslets^2; got $n_spots for n_lenslets=$n_sub"))
    gap_px = Int(gap)
    gap_px >= 0 || throw(ArgumentError("gap must be non-negative"))
    size(image) == (n_sub * n_y + (n_sub - 1) * gap_px, n_sub * n_x + (n_sub - 1) * gap_px) ||
        throw(DimensionMismatch("output image has size $(size(image)); expected " *
                                 "$((n_sub * n_y + (n_sub - 1) * gap_px, n_sub * n_x + (n_sub - 1) * gap_px))"))

    _shack_hartmann_detector_image!(image_style, spot_style, image, spot_cube, n_sub, n_y, n_x, gap_px, gap_value)
    return image
end

function _shack_hartmann_detector_image!(::ExecutionStyle, ::ExecutionStyle, image::AbstractMatrix,
    spot_cube::AbstractArray, n_sub::Int, n_y::Int, n_x::Int, gap::Int, gap_value)
    throw(InvalidConfiguration("Shack-Hartmann detector image output and spot cube must use the same execution backend"))
end

function _shack_hartmann_detector_image!(::ScalarCPUStyle, ::ScalarCPUStyle, image::AbstractMatrix,
    spot_cube::AbstractArray, n_sub::Int, n_y::Int, n_x::Int, gap::Int, gap_value)
    T = eltype(image)
    fill!(image, detector_output_value(T, gap_value))
    @inbounds for i in 1:n_sub, j in 1:n_sub
        idx = (i - 1) * n_sub + j
        y0 = (i - 1) * (n_y + gap) + 1
        x0 = (j - 1) * (n_x + gap) + 1
        for jj in 1:n_x, ii in 1:n_y
            image[y0 + ii - 1, x0 + jj - 1] = detector_output_value(T, spot_cube[idx, ii, jj])
        end
    end
    return image
end

function _shack_hartmann_detector_image!(style::AcceleratorStyle{B}, ::AcceleratorStyle{B},
    image::AbstractMatrix, spot_cube::AbstractArray, n_sub::Int, n_y::Int, n_x::Int, gap::Int, gap_value) where {B}
    launch_kernel!(style, shack_hartmann_detector_image_kernel!, image, spot_cube, n_sub, n_y, n_x, gap,
        eltype(image)(gap_value); ndrange=size(image))
    return image
end

@inline supports_valid_subaperture_mask(wfs::AbstractWFS) = !isnothing(valid_subaperture_mask(wfs))
@inline supports_reference_signal(wfs::AbstractWFS) = !isnothing(reference_signal(wfs))
@inline supports_camera_frame(wfs::AbstractWFS) = !isnothing(camera_frame(wfs))

@inline valid_subaperture_mask(wfs::ShackHartmannWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::ShackHartmannWFS) = wfs.state.reference_signal_2d
@inline shack_hartmann_detector_image(wfs::ShackHartmannWFS; kwargs...) =
    shack_hartmann_detector_image(sh_exported_spot_cube(wfs), wfs.params.n_lenslets; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, ::Nothing; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, det::AbstractDetector; kwargs...) =
    shack_hartmann_detector_image(wfs; output_type=detector_output_type(det), kwargs...)

@inline valid_subaperture_mask(wfs::PyramidWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::PyramidWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::PyramidWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::BioEdgeWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::BioEdgeWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::BioEdgeWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::ZernikeWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::ZernikeWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::ZernikeWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::CurvatureWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::CurvatureWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::CurvatureWFS) = wfs.state.camera_frame

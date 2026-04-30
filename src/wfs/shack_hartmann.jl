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

@inline valid_subaperture_mask(wfs::ShackHartmannWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::ShackHartmannWFS) = wfs.state.reference_signal_2d

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

@inline shack_hartmann_detector_image(wfs::ShackHartmannWFS; kwargs...) =
    shack_hartmann_detector_image(sh_exported_spot_cube(wfs), wfs.params.n_lenslets; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, ::Nothing; kwargs...) = shack_hartmann_detector_image(wfs; kwargs...)
@inline wfs_detector_image(wfs::ShackHartmannWFS, det::AbstractDetector; kwargs...) =
    shack_hartmann_detector_image(wfs; output_type=detector_output_type(det), kwargs...)

@inline wfs_output_frame(wfs::ShackHartmannWFS{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame(wfs::ShackHartmannWFS{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmannWFS{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmannWFS{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline wfs_output_metadata(wfs::ShackHartmannWFS) = (
    n_lenslets=subaperture_layout(wfs).n_subap,
    n_valid_subap=n_valid_subapertures(subaperture_layout(wfs)),
    subap_pixels=subaperture_layout(wfs).subap_pixels,
    pitch_m=subaperture_layout(wfs).pitch_m,
    slopes_units=subaperture_calibration(wfs).slopes_units,
    calibrated=subaperture_calibration(wfs).calibrated,
)

@inline supports_prepared_runtime(::ShackHartmannWFS{<:Diffractive}, ::AbstractSource) = true
@inline supports_prepared_runtime(::ShackHartmannWFS{<:Diffractive}, ::Asterism) = true
@inline supports_detector_output(::ShackHartmannWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::ShackHartmannWFS, ::Asterism) = true
@inline supports_stacked_sources(::ShackHartmannWFS, ::SpectralSource) = true
@inline supports_stacked_sources(::ShackHartmannWFS, ::ExtendedSource) = true
@inline supports_grouped_execution(::ShackHartmannWFS{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(::ShackHartmannWFS{<:Diffractive}, ::SpectralSource) = true
@inline supports_grouped_execution(::ShackHartmannWFS{<:Diffractive}, ::ExtendedSource) = true

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, tel::Telescope, src::AbstractSource)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, tel::Telescope, src::SpectralSource)
    prepare_sampling!(wfs, tel, spectral_reference_source(src))
    ensure_sh_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmannWFS{<:Diffractive}, tel::Telescope, ast::Asterism)
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    return wfs
end

@inline function _measure_for_calibration!(wfs::ShackHartmannWFS{<:Diffractive}, tel::Telescope, src::AbstractSource)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    return measure!(wfs, tel, src)
end

@inline function _measure_for_calibration!(wfs::ShackHartmannWFS{<:Diffractive}, tel::Telescope, src::SpectralSource)
    prepare_sampling!(wfs, tel, spectral_reference_source(src))
    ensure_sh_calibration!(wfs, tel, src)
    return measure!(wfs, tel, src)
end

@inline function set_runtime_wfs_output_policy!(wfs::ShackHartmannWFS{<:Diffractive}, outputs)
    wfs.state.export_pixels_enabled = outputs.wfs_pixels
    return wfs
end

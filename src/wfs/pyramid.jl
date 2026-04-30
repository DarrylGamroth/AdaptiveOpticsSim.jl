#
# Pyramid wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("pyramid/setup.jl")
include("pyramid/measure.jl")
include("pyramid/optics.jl")
include("pyramid/signals.jl")

@inline valid_subaperture_mask(wfs::PyramidWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::PyramidWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::PyramidWFS) = wfs.state.camera_frame

@inline wfs_output_frame_prototype(wfs::PyramidWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::PyramidWFS, det::AbstractDetector) = camera_frame(wfs)

@inline supports_prepared_runtime(::PyramidWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::PyramidWFS, ::Asterism) = true
@inline supports_detector_output(::PyramidWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::PyramidWFS, ::Asterism) = true
@inline supports_stacked_sources(::PyramidWFS, ::SpectralSource) = true
@inline supports_stacked_sources(::PyramidWFS, ::ExtendedSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::SpectralSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::ExtendedSource) = true

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, ast.sources[1])
    return wfs
end

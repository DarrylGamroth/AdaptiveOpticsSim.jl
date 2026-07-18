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
include("pyramid/stages.jl")

@inline valid_subaperture_mask(wfs::PyramidWFS) = wfs.estimator.state.valid_mask
@inline reference_signal(wfs::PyramidWFS) = wfs.estimator.state.reference_signal_2d
@inline slopes(wfs::PyramidWFS) = wfs.estimator.state.slopes
@inline camera_frame(wfs::PyramidWFS{<:Diffractive}) =
    wfs.acquisition.state.camera_frame
@inline camera_frame(::PyramidWFS{<:Geometric}) = nothing

@inline wfs_output_frame_prototype(wfs::PyramidWFS{<:Diffractive},
    ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::PyramidWFS{<:Diffractive},
    det::AbstractDetector) = camera_frame(wfs)

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
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, common_source)
    return wfs
end

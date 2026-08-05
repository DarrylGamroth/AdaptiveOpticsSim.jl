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
@inline wfs_calibration_signature(wfs::PyramidWFS) =
    wfs.estimator.state.calibration_signature

@inline supports_prepared_runtime(::PyramidWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::PyramidWFS, ::Asterism) = true
@inline supports_detector_output(::PyramidWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::PyramidWFS, ::Asterism) = true
@inline supports_stacked_sources(::PyramidWFS, ::SpectralSource) = true
@inline supports_stacked_sources(::PyramidWFS, ::ExtendedSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::SpectralSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::ExtendedSource) = true

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    prepare_pyramid_sampling!(wfs, pupil)
    ensure_pyramid_calibration!(wfs, pupil, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource)
    prepare_pyramid_sampling!(wfs, pupil)
    ensure_pyramid_calibration!(wfs, pupil, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    prepare_pyramid_sampling!(wfs, pupil)
    ensure_pyramid_calibration!(wfs, pupil, common_source)
    return wfs
end

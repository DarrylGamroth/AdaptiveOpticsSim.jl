#
# Pyramid wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("pyramid/setup.jl")
include("pyramid/optics.jl")
include("pyramid/stages.jl")

@inline supports_prepared_runtime(::PyramidWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::PyramidWFS, ::Asterism) = true
@inline supports_detector_output(::PyramidWFS, ::AbstractDetector) = true
@inline supports_stacked_sources(::PyramidWFS, ::Asterism) = true
@inline supports_stacked_sources(::PyramidWFS, ::SpectralSource) = true
@inline supports_stacked_sources(::PyramidWFS, ::ExtendedSource) = true
@inline supports_grouped_execution(::PyramidWFS, ::Asterism) = true
@inline supports_grouped_execution(::PyramidWFS, ::SpectralSource) = true
@inline supports_grouped_execution(::PyramidWFS, ::ExtendedSource) = true

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    prepare_pyramid_sampling!(wfs, pupil)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource)
    prepare_pyramid_sampling!(wfs, pupil)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    prepare_pyramid_sampling!(wfs, pupil)
    return wfs
end

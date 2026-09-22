#
# Bi-O-edge wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("bi_o_edge/setup.jl")
include("bi_o_edge/optics.jl")
include("bi_o_edge/stages.jl")

@inline supports_prepared_runtime(::BiOEdgeWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::BiOEdgeWFS, ::Asterism) = true
@inline supports_detector_output(::BiOEdgeWFS, ::AbstractDetector) = true
@inline supports_stacked_sources(::BiOEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::BiOEdgeWFS, ::Asterism) = true

@inline function prepare_runtime_wfs!(wfs::BiOEdgeWFS,
    pupil::PupilFunction, ::AbstractSource)
    prepare_bi_o_edge_sampling!(wfs, pupil)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::BiOEdgeWFS,
    pupil::PupilFunction, ::Asterism)
    prepare_bi_o_edge_sampling!(wfs, pupil)
    return wfs
end

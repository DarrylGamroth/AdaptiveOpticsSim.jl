#
# Bi-O-edge wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("bi_o_edge/setup.jl")
include("bi_o_edge/measure.jl")
include("bi_o_edge/signals.jl")
include("bi_o_edge/stages.jl")

@inline valid_subaperture_mask(wfs::BiOEdgeWFS) = wfs.estimator.state.valid_mask
@inline reference_signal(wfs::BiOEdgeWFS) = wfs.estimator.state.reference_signal_2d
@inline slopes(wfs::BiOEdgeWFS) = bi_o_edge_estimator_products(wfs).slopes
@inline wfs_calibration_signature(wfs::BiOEdgeWFS) =
    wfs.estimator.state.calibration_signature

@inline supports_detector_output(::BiOEdgeWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::BiOEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::BiOEdgeWFS{<:Diffractive}, ::Asterism) = true

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
@inline slopes(wfs::BiOEdgeWFS) = wfs.estimator.state.slopes
@inline camera_frame(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.acquisition.state.camera_frame
@inline camera_frame(::BiOEdgeWFS{<:Geometric}) = nothing
@inline wfs_calibration_signature(wfs::BiOEdgeWFS) =
    wfs.estimator.state.calibration_signature

@inline wfs_output_frame_prototype(wfs::BiOEdgeWFS{<:Diffractive},
    ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::BiOEdgeWFS{<:Diffractive},
    det::AbstractDetector) = camera_frame(wfs)

@inline supports_detector_output(::BiOEdgeWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::BiOEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::BiOEdgeWFS{<:Diffractive}, ::Asterism) = true

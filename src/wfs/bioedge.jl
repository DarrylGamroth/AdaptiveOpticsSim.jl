#
# BioEdge wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("bioedge/setup.jl")
include("bioedge/measure.jl")
include("bioedge/signals.jl")

@inline valid_subaperture_mask(wfs::BioEdgeWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::BioEdgeWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::BioEdgeWFS) = wfs.state.camera_frame

@inline wfs_output_frame_prototype(wfs::BioEdgeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::BioEdgeWFS, det::AbstractDetector) = camera_frame(wfs)

@inline supports_detector_output(::BioEdgeWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::BioEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::BioEdgeWFS{<:Diffractive}, ::Asterism) = true

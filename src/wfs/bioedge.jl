#
# BioEdge wavefront sensing
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("bioedge/setup.jl")
include("bioedge/measure.jl")
include("bioedge/signals.jl")
include("bioedge/stages.jl")

@inline valid_subaperture_mask(wfs::BioEdgeWFS) = wfs.estimator.state.valid_mask
@inline reference_signal(wfs::BioEdgeWFS) = wfs.estimator.state.reference_signal_2d
@inline slopes(wfs::BioEdgeWFS) = wfs.estimator.state.slopes
@inline camera_frame(wfs::BioEdgeWFS{<:Diffractive}) =
    wfs.acquisition.state.camera_frame
@inline camera_frame(::BioEdgeWFS{<:Geometric}) = nothing

@inline wfs_output_frame_prototype(wfs::BioEdgeWFS{<:Diffractive},
    ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::BioEdgeWFS{<:Diffractive},
    det::AbstractDetector) = camera_frame(wfs)

@inline supports_detector_output(::BioEdgeWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_stacked_sources(::BioEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::BioEdgeWFS{<:Diffractive}, ::Asterism) = true

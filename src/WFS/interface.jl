#
# Wavefront-sensor interface accessors
#
# This file defines the maintained subsystem-level accessors for `AbstractWFS`.
# The goal is to keep runtime, calibration, and tests off direct `state` field
# access when the accessed quantity is part of the supported WFS contract.
#

"""
    slopes(wfs::AbstractWFS)

Return the current exported 1-D signal vector for a wavefront sensor.

For geometric and centroid-based sensors this is the slope vector. For sensors
such as `ZernikeWFS` and `CurvatureWFS`, the maintained runtime contract still
uses the `slopes` name even though the stored values represent normalized
signal samples rather than geometric slopes.
"""
@inline slopes(wfs::AbstractWFS) = wfs.state.slopes

"""
    valid_subaperture_mask(wfs::AbstractWFS)

Return the maintained valid-subaperture/support mask when the WFS exposes one.
Returns `nothing` for sensor families that do not export a stable valid-mask
surface.
"""
@inline valid_subaperture_mask(::AbstractWFS) = nothing

"""
    reference_signal(wfs::AbstractWFS)

Return the stored reference signal surface used by the WFS calibration path, or
`nothing` when the family does not expose one as part of the maintained
interface.
"""
@inline reference_signal(::AbstractWFS) = nothing

"""
    camera_frame(wfs::AbstractWFS)

Return the maintained internal camera/readout frame for detector-like WFS
families, or `nothing` when the family does not expose a stable camera-frame
surface.
"""
@inline camera_frame(::AbstractWFS) = nothing

@inline supports_valid_subaperture_mask(wfs::AbstractWFS) = !isnothing(valid_subaperture_mask(wfs))
@inline supports_reference_signal(wfs::AbstractWFS) = !isnothing(reference_signal(wfs))
@inline supports_camera_frame(wfs::AbstractWFS) = !isnothing(camera_frame(wfs))

@inline valid_subaperture_mask(wfs::ShackHartmann) = wfs.state.valid_mask
@inline reference_signal(wfs::ShackHartmann) = wfs.state.reference_signal_2d

@inline valid_subaperture_mask(wfs::PyramidWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::PyramidWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::PyramidWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::BioEdgeWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::BioEdgeWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::BioEdgeWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::ZernikeWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::ZernikeWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::ZernikeWFS) = wfs.state.camera_frame

@inline valid_subaperture_mask(wfs::CurvatureWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::CurvatureWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::CurvatureWFS) = wfs.state.camera_frame

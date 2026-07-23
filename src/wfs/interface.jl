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
function slopes end

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
    wfs_calibration_signature(wfs::AbstractWFS)

Return the optical calibration signature currently bound to a wavefront
sensor. This accessor intentionally hides the family-specific ownership of
calibration state.
"""
@inline wfs_calibration_signature(wfs::ShackHartmannWFS) =
    wfs.calibration.signature
@inline wfs_calibration_signature(wfs::PyramidWFS) =
    wfs.estimator.state.calibration_signature
@inline wfs_calibration_signature(wfs::BioEdgeWFS) =
    wfs.estimator.state.calibration_signature
@inline wfs_calibration_signature(wfs::ZernikeWFS) =
    wfs.estimator.state.calibration_signature
@inline wfs_calibration_signature(wfs::CurvatureWFS) =
    wfs.estimator.state.calibration_signature

"""
    camera_frame(wfs::AbstractWFS)

Return the maintained internal camera/readout frame for detector-like WFS
families, or `nothing` when the family does not expose a stable camera-frame
surface.
"""
@inline camera_frame(::AbstractWFS) = nothing

"""
    wfs_detector_image(wfs)
    wfs_detector_image(wfs, det)

Return the detector-image product for a wavefront sensor.

For WFS families with a native 2-D camera frame, this returns `camera_frame(wfs)`
when no detector is provided and `output_frame(det)` after detector-coupled
measurement when a detector is provided. Shack-Hartmann WFS objects store
detector-coupled lenslet spots as a cube, so this accessor returns a tiled 2-D
detector mosaic.
"""
function wfs_detector_image(wfs::AbstractWFS)
    frame = camera_frame(wfs)
    isnothing(frame) && throw(InvalidConfiguration("WFS of type $(typeof(wfs)) does not expose a detector image"))
    return frame
end

@inline wfs_detector_image(::AbstractWFS, det::AbstractDetector) = output_frame(det)

@inline supports_valid_subaperture_mask(wfs::AbstractWFS) = !isnothing(valid_subaperture_mask(wfs))
@inline supports_reference_signal(wfs::AbstractWFS) = !isnothing(reference_signal(wfs))
@inline supports_camera_frame(wfs::AbstractWFS) = !isnothing(camera_frame(wfs))

@inline wfs_output_frame(wfs::AbstractWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame(::AbstractWFS, detector::AbstractDetector) =
    output_frame(detector)
@inline wfs_output_frame_prototype(wfs::AbstractWFS, detector) =
    wfs_output_frame(wfs, detector)
@inline wfs_output_metadata(::AbstractWFS) = nothing

"""
    supports_prepared_runtime(wfs, source)

Return whether a WFS/source pairing provides meaningful preparation for
repeated measurement.
"""
@inline supports_prepared_runtime(::AbstractWFS, ::Any) = false

"""Return whether `wfs` supports detector-coupled measurement with `detector`."""
@inline supports_detector_output(::AbstractWFS, ::AbstractDetector) = false

"""Return whether `wfs` supports a source represented by stacked components."""
@inline supports_stacked_sources(::AbstractWFS, ::Any) = false

"""Return whether `wfs` provides grouped execution for `source`."""
@inline supports_grouped_execution(::AbstractWFS, ::Any) = false

"""
    prepare_runtime_wfs!(wfs, pupil, source)

Perform WFS-specific preparation for repeated measurement. The default is a
no-op. This is a WFS preparation hook, not a simulation-runtime constructor.
"""
@inline prepare_runtime_wfs!(wfs::AbstractWFS, ::PupilFunction, source) = wfs

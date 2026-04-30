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

@inline function prepropagate_runtime_wfs!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource, rng::AbstractRNG)
    advance!(atm, tel, rng)
    propagate!(atm, tel)
    apply!(optic, tel, DMAdditive())
    return nothing
end

@inline function measure_runtime_wfs!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource, rng::AbstractRNG)
    measure!(wfs, tel, src)
    return nothing
end

@inline function measure_runtime_wfs!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    measure!(wfs, tel, src, det; rng=rng)
    return nothing
end

@inline function finish_runtime_wfs_sensing!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource)
    return nothing
end

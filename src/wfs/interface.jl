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

For geometric and centroid-based sensors this is the slope vector. Other
families may use the maintained `slopes` contract for normalized signal
samples that are not geometric slopes.
"""
function slopes end

"""
    measure!(wfs, pupil[, source][, detector]; rng)

Update `wfs` from a pupil-plane input and return its maintained measurement
product. Concrete wavefront-sensor families define the supported source and
detector combinations.
"""
function measure! end

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
function wfs_calibration_signature end

"""
    wfs_detector_image(wfs, det)

Return the current two-dimensional output storage of a detector compatible
with `wfs`. The caller must first perform detector-coupled acquisition. Use
`Detectors.output_frame(det)` directly for non-image detector products.
"""
function wfs_detector_image(wfs::AbstractWFS, det::AbstractDetector)
    supports_detector_output(wfs, det) ||
        throw(InvalidConfiguration(
            "WFS of type $(typeof(wfs)) does not support detector output " *
            "from $(typeof(det))"))
    frame = output_frame(det)
    ndims(frame) == 2 ||
        throw(InvalidConfiguration(
            "detector output for WFS of type $(typeof(wfs)) is not a " *
            "two-dimensional image; use Detectors.output_frame(det)"))
    return frame
end

@inline supports_valid_subaperture_mask(wfs::AbstractWFS) = !isnothing(valid_subaperture_mask(wfs))
@inline supports_reference_signal(wfs::AbstractWFS) = !isnothing(reference_signal(wfs))
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

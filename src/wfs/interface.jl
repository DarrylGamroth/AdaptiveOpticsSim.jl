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
    wfs.state.calibration_signature
@inline wfs_calibration_signature(wfs::CurvatureWFS) =
    wfs.state.calibration_signature

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

@inline propagate_runtime_atmosphere!(atm::AbstractAtmosphere, tel::Telescope,
    src::AbstractSource) = propagate!(atm, tel, src)
@inline propagate_runtime_atmosphere!(atm::AbstractAtmosphere, tel::Telescope,
    ::Asterism) = propagate!(atm, tel)

@inline prepare_runtime_atmosphere_path(::AbstractAtmosphere,
    ::Telescope, ::AbstractSource) = nothing

@inline prepare_runtime_atmosphere_path(atm::AbstractTimedAtmosphere,
    tel::Telescope, src::AbstractSource) =
    prepare_atmosphere_renderer(atm, tel, src)

@inline render_prepared_atmosphere_path!(::Nothing,
    atm::AbstractAtmosphere, tel::Telescope, src::AbstractSource) =
    propagate_runtime_atmosphere!(atm, tel, src)

@inline function render_prepared_atmosphere_path!(
    renderer::AtmosphereDirectionRenderer,
    atm::AbstractTimedAtmosphere, tel::Telescope, ::AbstractSource)
    render_atmosphere_opd!(opd_map(tel), renderer, atm,
        current_epoch(atm))
    return tel
end

@inline function advance_runtime_atmosphere!(wfs::AbstractWFS,
    atm::AbstractAtmosphere, tel::Telescope, atmosphere_step::Real,
    rng::AbstractRNG)
    advance!(atm, tel, rng)
    return nothing
end

@inline function advance_runtime_atmosphere!(::AbstractWFS,
    atm::AbstractTimedAtmosphere, ::Telescope, atmosphere_step::Real,
    rng::AbstractRNG)
    advance_by!(atm, atmosphere_step, rng)
    return nothing
end

@inline function render_runtime_wfs_path!(wfs::AbstractWFS,
    atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    src::AbstractSource)
    propagate_runtime_atmosphere!(atm, tel, src)
    apply!(optic, tel, DMAdditive())
    return nothing
end

@inline function render_runtime_wfs_path!(renderer, wfs::AbstractWFS,
    atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, src::AbstractSource)
    render_prepared_atmosphere_path!(renderer, atm, tel, src)
    apply!(optic, tel, DMAdditive())
    return nothing
end

@inline function prepropagate_runtime_wfs!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource,
    atmosphere_step::Real, rng::AbstractRNG)
    advance_runtime_atmosphere!(wfs, atm, tel, atmosphere_step, rng)
    render_runtime_wfs_path!(wfs, atm, tel, optic, src)
    return nothing
end

@inline prepare_shared_runtime_wfs!(wfs::AbstractWFS, atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource) = nothing

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
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource)
    return nothing
end

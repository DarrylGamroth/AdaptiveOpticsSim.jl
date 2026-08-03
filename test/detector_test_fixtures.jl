function detector_test_intensity_map(values::AbstractMatrix{T};
    kind::AbstractOpticalPlaneKind=FocalPlane(),
    sampling::NTuple{2,T}=(one(T), one(T)),
    normalization::AbstractOpticalNormalization=PhotonRateNormalization(),
    spatial_measure::AbstractSpatialMeasure=CellIntegratedMeasure(),
    coherence::AbstractCombinationPolicy=IncoherentIntensityAddition(),
    spectral::AbstractSpectralCoordinate=MonochromaticChannel(0.55e-6)) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(kind, values;
        coordinate_domain=AngularCoordinates(), sampling=sampling,
        normalization=normalization, spatial_measure=spatial_measure,
        coherence=coherence, spectral=spectral)
    return IntensityMap(metadata, values)
end

function detector_state_snapshot(det::Detector)
    names = fieldnames(typeof(det.state))
    values = map(name -> getfield(det.state, name), names)
    return NamedTuple{names}(values)
end

function detector_state_matches_snapshot(det::Detector, snapshot::NamedTuple)
    return all(name -> getfield(det.state, name) === getfield(snapshot, name),
        keys(snapshot))
end

struct UnkeyedCalibrationFrameResponse{T<:AbstractFloat} <:
    AdaptiveOpticsSim.Detectors.AbstractFrameResponse
    alpha::T
end

struct UnsupportedExactTargetReadoutProducts <:
    AdaptiveOpticsSim.Detectors.FrameReadoutProducts
end

function AdaptiveOpticsSim.Detectors.convert_frame_response_model(
    model::UnkeyedCalibrationFrameResponse, ::Type{T}, backend) where
    {T<:AbstractFloat}
    return UnkeyedCalibrationFrameResponse{T}(T(model.alpha))
end

AdaptiveOpticsSim.Detectors.validate_frame_response_model(
    model::UnkeyedCalibrationFrameResponse) = model

macro test_detector_allocation(expression)
    return quote
        if coverage_instrumented()
            @test_skip "detector allocation assertion is disabled under coverage instrumentation"
        else
            @test $(esc(expression))
        end
    end
end

function prepared_detector_capture_allocations(prepared, rng)
    capture!(prepared, rng)
    return @allocated capture!(prepared, rng)
end

function prepared_detector_readiness_allocations(prepared)
    AdaptiveOpticsSim.Detectors._require_prepared_whole_acquisition(prepared)
    return @allocated AdaptiveOpticsSim.Detectors._require_prepared_whole_acquisition(
        prepared)
end

function prepared_first_detector_capture_allocations(builder, map)
    warm_detector = builder()
    warm_plan = prepare_detector_acquisition(warm_detector, map)
    capture!(warm_plan, Xoshiro(2400))

    detector = builder()
    plan = prepare_detector_acquisition(detector, map)
    rng = Xoshiro(2401)
    return @allocated capture!(plan, rng)
end

function prepared_detector_exposed_storage_is_zero(det::Detector)
    products = (
        output_frame(det),
        detector_reference_frame(det),
        detector_signal_frame(det),
        detector_combined_frame(det),
        detector_reference_cube(det),
        detector_signal_cube(det),
        detector_read_cube(det),
        detector_read_times(det),
        detector_ramp_slope(det),
        detector_ramp_intercept(det),
    )
    return all(product -> isnothing(product) || all(iszero, product),
        products)
end

function prepared_incremental_capture_allocations(prepared, rng,
    integration_duration)
    capture!(prepared; rng=rng, integration_duration=integration_duration)
    return @allocated capture!(prepared; rng=rng,
        integration_duration=integration_duration)
end

function prepared_first_incremental_capture_allocations(det, prepared, rng,
    integration_duration)
    capture!(prepared; rng=rng, integration_duration=integration_duration)
    reset_integration!(det)
    return @allocated capture!(prepared; rng=rng,
        integration_duration=integration_duration)
end

function fixed_stack_capture_allocations(det, cube, scratch, rng)
    AdaptiveOpticsSim.Detectors.capture_stack!(det, cube, scratch, rng)
    fill!(cube, one(eltype(cube)))
    return @allocated AdaptiveOpticsSim.Detectors.capture_stack!(det, cube, scratch,
        rng)
end

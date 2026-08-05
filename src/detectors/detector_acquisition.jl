"""
    DetectorAcquisitionPlan

Run-immutable numerical contract for frame-detector acquisition. The plan owns
validated configuration, optical-plane description, dimensions, and
radiometric scaling. It does not own or identify live detector state, scratch,
products, or mutable intensity storage. Construct it through
[`prepare_detector_acquisition`](@ref), which returns the exact prepared owner
used for execution.
"""
struct _DetectorAcquisitionPlanToken end
const _DETECTOR_ACQUISITION_PLAN_TOKEN = _DetectorAcquisitionPlanToken()

struct DetectorAcquisitionPlan{
    P<:DetectorParams,
    M<:OpticalPlaneMetadata,
    T<:AbstractFloat,
}
    detector_params::P
    input_metadata::M
    input_shape::Tuple{Int,Int}
    frame_shape::Tuple{Int,Int}
    output_shape::Tuple{Int,Int}
    rate_scale::T
    quantum_efficiency::T

    function DetectorAcquisitionPlan(::_DetectorAcquisitionPlanToken,
        detector_params::P, input_metadata::M,
        input_shape::Tuple{Int,Int}, frame_shape::Tuple{Int,Int},
        output_shape::Tuple{Int,Int}, rate_scale::T,
        quantum_efficiency::T) where {
        P<:DetectorParams,M<:OpticalPlaneMetadata,T<:AbstractFloat,
    }
        return new{P,M,T}(detector_params, input_metadata, input_shape,
            frame_shape, output_shape, rate_scale, quantum_efficiency)
    end
end

struct _PreparedDetectorAcquisitionToken end
const _PREPARED_DETECTOR_ACQUISITION_TOKEN =
    _PreparedDetectorAcquisitionToken()

"""
    PreparedDetectorAcquisition

Exact single-writer owner for repeated frame-detector acquisition. It binds one
detector, one intensity product, the run-immutable plan, persistent detector
state, replaceable workspace, caller-visible products, and snapshots of every
replaceable owner member. Construct it with
[`prepare_detector_acquisition`](@ref) and execute it with `capture!`.
"""
struct PreparedDetectorAcquisition{
    D<:Detector,
    I<:IntensityMap,
    P<:DetectorAcquisitionPlan,
    S<:DetectorState,
    W<:DetectorWorkspace,
    O<:DetectorProducts,
    SB<:NamedTuple,
    WB<:NamedTuple,
    PB<:NamedTuple,
}
    detector::D
    input::I
    plan::P
    state::S
    workspace::W
    products::O
    state_binding::SB
    workspace_binding::WB
    product_binding::PB

    function PreparedDetectorAcquisition(
        ::_PreparedDetectorAcquisitionToken, detector::D, input::I,
        plan::P, state::S, workspace::W, products::O,
        state_binding::SB, workspace_binding::WB,
        product_binding::PB) where {
        D<:Detector,I<:IntensityMap,P<:DetectorAcquisitionPlan,
        S<:DetectorState,W<:DetectorWorkspace,O<:DetectorProducts,
        SB<:NamedTuple,WB<:NamedTuple,PB<:NamedTuple,
    }
        return new{D,I,P,S,W,O,SB,WB,PB}(detector, input, plan, state,
            workspace, products, state_binding, workspace_binding,
            product_binding)
    end
end

@inline detector_acquisition_detector(prepared::PreparedDetectorAcquisition) =
    prepared.detector
@inline detector_acquisition_input(prepared::PreparedDetectorAcquisition) =
    prepared.input
@inline detector_acquisition_plan(prepared::PreparedDetectorAcquisition) =
    prepared.plan
@inline detector_acquisition_state(prepared::PreparedDetectorAcquisition) =
    prepared.state
@inline detector_acquisition_workspace(
    prepared::PreparedDetectorAcquisition) = prepared.workspace
@inline detector_acquisition_products(prepared::PreparedDetectorAcquisition) =
    prepared.products

@inline _require_detector_acquisition_plane(::FocalPlane) = nothing
@inline _require_detector_acquisition_plane(::DetectorPlane) = nothing

function _require_detector_acquisition_plane(::AbstractOpticalPlaneKind)
    throw(InvalidConfiguration(
        "detector acquisition requires a focal-plane or detector-plane intensity map"))
end

@inline function _normalization_rate_scale(::PhotonRateNormalization,
    ::Nothing, ::Type{T}) where {T<:AbstractFloat}
    return one(T)
end

function _normalization_rate_scale(::PhotonRateNormalization,
    normalized_to_photon_rate::Real, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "normalized_to_photon_rate is only valid for dimensionless intensity maps"))
end

function _normalization_rate_scale(::DimensionlessNormalization,
    ::Nothing, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "dimensionless intensity requires an explicit normalized_to_photon_rate scale"))
end

function _normalization_rate_scale(::DimensionlessNormalization,
    normalized_to_photon_rate::Real, ::Type{T}) where {T<:AbstractFloat}
    scale = T(normalized_to_photon_rate)
    isfinite(scale) && scale > zero(T) || throw(InvalidConfiguration(
        "normalized_to_photon_rate must be finite and > 0"))
    return scale
end

function _normalization_rate_scale(::AbstractOpticalNormalization,
    normalized_to_photon_rate::Union{Nothing,Real},
    ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires photon-rate or explicitly scaled dimensionless normalization"))
end

@inline _spatial_rate_scale(::CellIntegratedMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat} = one(T)

function _spatial_rate_scale(::SpatialDensityMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat}
    dx = T(metadata.sampling[1])
    dy = T(metadata.sampling[2])
    isfinite(dx) && dx > zero(T) && isfinite(dy) && dy > zero(T) ||
        throw(InvalidConfiguration(
            "spatial-density acquisition requires finite positive plane sampling"))
    cell_measure = dx * dy
    isfinite(cell_measure) && cell_measure > zero(T) ||
        throw(InvalidConfiguration(
            "spatial-density plane sampling has an unrepresentable cell measure"))
    return cell_measure
end

function _spatial_rate_scale(::AbstractSpatialMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires cell-integrated or spatial-density intensity samples"))
end

@inline _require_detector_incoherent(::IncoherentIntensityAddition) = nothing

function _require_detector_incoherent(::AbstractCombinationPolicy)
    throw(InvalidConfiguration(
        "detector intensity acquisition requires incoherent-intensity semantics"))
end

@inline function _prepared_quantum_efficiency(det::Detector,
    channel::MonochromaticChannel, ::Type{T}) where {T<:AbstractFloat}
    return T(qe_at(det.params.quantum_efficiency_model,
        channel.wavelength_m))
end

@inline function _prepared_quantum_efficiency(det::Detector,
    ::IntegratedSpectralChannel, ::Type{T}) where {T<:AbstractFloat}
    return _integrated_spectral_qe(det.params.quantum_efficiency_model, T)
end

@inline _integrated_spectral_qe(model::ScalarQuantumEfficiency,
    ::Type{T}) where {T<:AbstractFloat} = T(model.value)

function _integrated_spectral_qe(::SampledQuantumEfficiency,
    ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "wavelength-dependent detector QE cannot be applied after passband integration"))
end

function _prepared_quantum_efficiency(det::Detector,
    ::AbstractSpectralCoordinate, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires a declared monochromatic or integrated spectral channel"))
end

@inline _require_prepared_response_sampling(::NullFrameResponse, ::Int) = nothing

function _require_prepared_response_sampling(::AbstractFrameResponse,
    psf_sampling::Int)
    psf_sampling == 1 || throw(InvalidConfiguration(
        "a non-null detector response currently requires psf_sampling == 1; " *
        "prepare an explicit optical-grid mapping before detector acquisition"))
    return nothing
end


function _require_finite_nonnegative_intensity(values::AbstractMatrix)
    isempty(values) && throw(InvalidConfiguration(
        "detector intensity input must not be empty"))
    # Preparation is allowed to synchronize. Validating through a host view
    # avoids backend-specific reduction compilation in this fallible setup
    # path and keeps repeated prepared acquisition device resident.
    minimum_value, maximum_value = extrema(host_array(values))
    isfinite(minimum_value) && isfinite(maximum_value) &&
        minimum_value >= zero(minimum_value) || throw(InvalidConfiguration(
            "detector intensity input values must be finite and nonnegative"))
    return nothing
end

@inline _require_detector_defect_shape(::NullDetectorDefectModel,
    ::Tuple{Int,Int}) = nothing

function _require_detector_defect_shape(model::PixelResponseNonuniformity,
    frame_shape::Tuple{Int,Int})
    size(model.gain_map) == frame_shape || throw(DimensionMismatchError(
        "PixelResponseNonuniformity gain_map size must match detector frame size"))
    return nothing
end

function _require_detector_defect_shape(model::DarkSignalNonuniformity,
    frame_shape::Tuple{Int,Int})
    size(model.dark_map) == frame_shape || throw(DimensionMismatchError(
        "DarkSignalNonuniformity dark_map size must match detector frame size"))
    return nothing
end

function _require_detector_defect_shape(model::BadPixelMask,
    frame_shape::Tuple{Int,Int})
    size(model.mask) == frame_shape || throw(DimensionMismatchError(
        "BadPixelMask mask size must match detector frame size"))
    return nothing
end

@inline _require_detector_defect_shapes(::Tuple{},
    ::Tuple{Int,Int}) = nothing

@inline function _require_detector_defect_shapes(models::Tuple,
    frame_shape::Tuple{Int,Int})
    _require_detector_defect_shape(first(models), frame_shape)
    _require_detector_defect_shapes(Base.tail(models), frame_shape)
    return nothing
end

@inline function _require_detector_defect_shape(
    model::CompositeDetectorDefectModel, frame_shape::Tuple{Int,Int})
    return _require_detector_defect_shapes(model.stages, frame_shape)
end

@inline _require_background_flux_shape(::NoBackground,
    ::Tuple{Int,Int}) = nothing
@inline _require_background_flux_shape(::ScalarBackground,
    ::Tuple{Int,Int}) = nothing

function _require_background_flux_shape(background::BackgroundFrame,
    frame_shape::Tuple{Int,Int})
    size(background.map) == frame_shape || throw(DimensionMismatchError(
        "background_flux size must match detector frame size"))
    return nothing
end

@inline _require_background_map_shape(::NoBackground,
    ::Tuple{Int,Int}) = nothing
@inline _require_background_map_shape(::ScalarBackground,
    ::Tuple{Int,Int}) = nothing

function _require_background_map_shape(background::BackgroundFrame,
    frame_shape::Tuple{Int,Int})
    size(background.map) == frame_shape || throw(DimensionMismatchError(
        "background_map size must match detector frame size"))
    return nothing
end

@inline _require_cmos_readout_shape(::NullCMOSReadNoise,
    ::Tuple{Int,Int}) = nothing

function _require_cmos_readout_shape(model::CMOSReadNoiseMap,
    frame_shape::Tuple{Int,Int})
    size(model.sigma) == frame_shape || throw(DimensionMismatchError(
        "CMOSReadNoiseMap sigma size must match detector frame size"))
    return nothing
end

@inline _require_cmos_output_shape(::NullCMOSOutputModel,
    ::Tuple{Int,Int}) = nothing

function _require_cmos_output_shape(model::StaticCMOSOutputPattern,
    frame_shape::Tuple{Int,Int})
    length(model.gains) * model.output_cols >= frame_shape[2] ||
        throw(DimensionMismatchError(
            "StaticCMOSOutputPattern does not cover detector columns"))
    return nothing
end

@inline _require_sensor_frame_shape(::AbstractFrameSensor,
    ::Tuple{Int,Int}) = nothing

@inline function _require_sensor_frame_shape(sensor::CMOSSensor,
    frame_shape::Tuple{Int,Int})
    _require_cmos_readout_shape(sensor.readout_noise_model, frame_shape)
    _require_cmos_output_shape(sensor.output_model, frame_shape)
    return nothing
end

function _require_detector_output_configuration(det::Detector)
    output = det.products.output_buffer
    output_host = det.workspace.output_buffer_host
    requires_output = det.params.output_type !== nothing ||
        det.params.readout_window !== nothing
    if requires_output
        output === nothing && throw(InvalidConfiguration(
            "Detector output type or readout window requires allocated output storage"))
        output_host === nothing && throw(InvalidConfiguration(
            "Detector output storage requires a prepared host output buffer"))
        return nothing
    end
    output === nothing || throw(InvalidConfiguration(
        "Detector output storage requires an output type or readout window"))
    output_host === nothing || throw(InvalidConfiguration(
        "Detector host output storage requires detector output storage"))
    return nothing
end

function _require_prepared_detector_storage(det::Detector,
    frame_shape::Tuple{Int,Int}, output_shape::Tuple{Int,Int})
    size(det.products.frame) == frame_shape || throw(DimensionMismatchError(
        "detector frame storage size must match prepared frame size"))
    output = det.products.output_buffer
    output === nothing && return nothing
    size(output) == output_shape || throw(DimensionMismatchError(
        "detector output storage size must match prepared readout size"))
    output_host = det.workspace.output_buffer_host
    output_host === nothing && throw(InvalidConfiguration(
        "Detector output storage requires a prepared host output buffer"))
    size(output_host) == output_shape || throw(DimensionMismatchError(
        "detector host output storage size must match prepared readout size"))
    return nothing
end

@inline _require_sensor_sampling_configuration(::AbstractFrameSensor,
    ::Tuple{Int,Int}, ::Union{Nothing,FrameWindow}, ::Real,
    ::Type{<:AbstractFloat}) = nothing

@inline function _require_sensor_sampling_configuration(
    sensor::AbstractHgCdTeSensor, frame_shape::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, exposure_duration::Real,
    ::Type{T}) where {T<:AbstractFloat}
    return _require_hgcdte_sampling_configuration(
        multi_read_sampling_mode(sensor),
        sensor, frame_shape, window, exposure_duration, T)
end

@inline _require_hgcdte_sampling_configuration(::FrameSamplingMode,
    ::AbstractHgCdTeSensor, ::Tuple{Int,Int},
    ::Union{Nothing,FrameWindow}, ::Real, ::Type{<:AbstractFloat}) = nothing

@inline function _require_hgcdte_sampling_configuration(
    mode::UpTheRampSampling, sensor::AbstractHgCdTeSensor,
    frame_shape::Tuple{Int,Int}, window::Union{Nothing,FrameWindow},
    exposure_duration::Real, ::Type{T}) where {T<:AbstractFloat}
    validate_up_the_ramp_schedule(sensor, frame_shape, window, mode,
        exposure_duration, T)
    return nothing
end

@inline function _require_detector_configuration(det::Detector,
    frame_shape::Tuple{Int,Int})
    _require_detector_defect_shape(det.params.defect_model, frame_shape)
    _require_background_flux_shape(det.background_flux, frame_shape)
    _require_background_map_shape(det.background_map, frame_shape)
    _require_sensor_frame_shape(det.params.sensor, frame_shape)
    _require_detector_output_configuration(det)
    _require_sensor_sampling_configuration(det.params.sensor, frame_shape,
        det.params.readout_window, det.params.exposure_duration,
        eltype(det.products.frame))
    return nothing
end

@inline _frame_readout_workspace_members(::NoFrameReadoutWorkspace) = ()
@inline _frame_readout_workspace_members(workspace::SkipperReadoutWorkspace) =
    (workspace.baseline_frame, workspace.sample_sum)
@inline function _frame_readout_workspace_members(
    workspace::MultiReadFrameReadoutWorkspace)
    return (workspace.reference_average, workspace.signal_average,
        workspace.reference_cube, workspace.signal_cube)
end
@inline function _frame_readout_workspace_members(
    workspace::UpTheRampReadoutWorkspace)
    return (workspace.slope, workspace.intercept, workspace.integrated,
        workspace.cube)
end

@inline _frame_readout_product_members(::NoFrameReadoutProducts) = ()
@inline _frame_readout_product_members(products::SkipperReadoutProducts) =
    (products.mean_frame, products.sample_count)
@inline function _frame_readout_product_members(
    products::SampledFrameReadoutProducts)
    return (products.reference_frame, products.signal_frame,
        products.read_cube)
end
@inline function _frame_readout_product_members(
    products::MultiReadFrameReadoutProducts)
    return (products.reference_frame, products.signal_frame,
        products.combined_frame, products.reference_cube,
        products.signal_cube, products.read_cube, products.read_offsets_s)
end
@inline function _frame_readout_product_members(
    products::UpTheRampReadoutProducts)
    return (products.slope_frame, products.intercept_frame,
        products.integrated_frame, products.read_cube, products.read_offsets_s,
        products.acquisition_kind)
end

@inline function _detector_state_members(state::DetectorState)
    return (state.accum_buffer, state.latent_buffer, state.thermal_state)
end

@inline function _detector_workspace_members(workspace::DetectorWorkspace)
    return (workspace.presampling_buffer, workspace.presampling_scratch,
        workspace.response_buffer, workspace.bin_buffer,
        workspace.temporal_buffer, workspace.noise_buffer,
        workspace.noise_buffer_host, workspace.batched_buffer_host,
        workspace.output_buffer_host, workspace.readout,
        _frame_readout_workspace_members(workspace.readout)...)
end

@inline function _detector_product_members(products::DetectorProducts)
    return (products.frame, products.output_buffer, products.readout,
        _frame_readout_product_members(products.readout)...)
end

@inline _storage_mightalias(value, other) = false

@inline function _storage_mightalias(value::AbstractArray,
    other::AbstractArray)
    return Base.mightalias(value, other)
end

@inline _storage_mightalias_any(value, ::Tuple{}) = false

@inline function _storage_mightalias_any(value, members::Tuple)
    return _storage_mightalias(value, first(members)) ||
        _storage_mightalias_any(value, Base.tail(members))
end

@inline function _detector_storage_mightalias(det::Detector,
    value::AbstractArray)
    return _storage_mightalias_any(value,
            _detector_state_members(det.state)) ||
        _storage_mightalias_any(value,
            _detector_workspace_members(det.workspace)) ||
        _storage_mightalias_any(value,
            _detector_product_members(det.products))
end

@inline _frame_readout_workspace_binding(::NoFrameReadoutWorkspace) = (;)
@inline _frame_readout_workspace_binding(workspace::SkipperReadoutWorkspace) =
    (baseline_frame=workspace.baseline_frame, sample_sum=workspace.sample_sum)
@inline function _frame_readout_workspace_binding(
    workspace::MultiReadFrameReadoutWorkspace)
    return (reference_average=workspace.reference_average,
        signal_average=workspace.signal_average,
        reference_cube=workspace.reference_cube,
        signal_cube=workspace.signal_cube)
end
@inline function _frame_readout_workspace_binding(
    workspace::UpTheRampReadoutWorkspace)
    return (slope=workspace.slope, intercept=workspace.intercept,
        integrated=workspace.integrated, cube=workspace.cube)
end

@inline _frame_readout_product_binding(::NoFrameReadoutProducts) = (;)
@inline _frame_readout_product_binding(products::SkipperReadoutProducts) =
    (mean_frame=products.mean_frame, sample_count=products.sample_count)
@inline function _frame_readout_product_binding(
    products::SampledFrameReadoutProducts)
    return (reference_frame=products.reference_frame,
        signal_frame=products.signal_frame, read_cube=products.read_cube)
end
@inline function _frame_readout_product_binding(
    products::MultiReadFrameReadoutProducts)
    return (reference_frame=products.reference_frame,
        signal_frame=products.signal_frame,
        combined_frame=products.combined_frame,
        reference_cube=products.reference_cube,
        signal_cube=products.signal_cube, read_cube=products.read_cube,
        read_offsets_s=products.read_offsets_s)
end
@inline function _frame_readout_product_binding(
    products::UpTheRampReadoutProducts)
    return (slope_frame=products.slope_frame,
        intercept_frame=products.intercept_frame,
        integrated_frame=products.integrated_frame,
        read_cube=products.read_cube, read_offsets_s=products.read_offsets_s)
end

@inline function _detector_state_binding(state::DetectorState)
    return (accum_buffer=state.accum_buffer,
        latent_buffer=state.latent_buffer, thermal_state=state.thermal_state)
end

@inline function _detector_workspace_binding(workspace::DetectorWorkspace)
    return (presampling_buffer=workspace.presampling_buffer,
        presampling_scratch=workspace.presampling_scratch,
        response_buffer=workspace.response_buffer,
        bin_buffer=workspace.bin_buffer,
        temporal_buffer=workspace.temporal_buffer,
        noise_buffer=workspace.noise_buffer,
        noise_buffer_host=workspace.noise_buffer_host,
        batched_buffer_host=workspace.batched_buffer_host,
        output_buffer_host=workspace.output_buffer_host,
        readout=workspace.readout,
        readout_binding=_frame_readout_workspace_binding(workspace.readout))
end

@inline function _detector_product_binding(products::DetectorProducts)
    return (frame=products.frame, output_buffer=products.output_buffer,
        readout=products.readout,
        readout_binding=_frame_readout_product_binding(products.readout))
end

@inline function _require_exact_binding(actual, expected,
    label::AbstractString)
    actual === expected || throw(InvalidConfiguration(
        "$label changed after detector acquisition preparation"))
    return nothing
end

@inline _require_frame_readout_workspace_binding(
    ::NoFrameReadoutWorkspace, ::NamedTuple{()}) = nothing

@inline function _require_frame_readout_workspace_binding(
    workspace::SkipperReadoutWorkspace, binding::NamedTuple)
    _require_exact_binding(workspace.baseline_frame, binding.baseline_frame,
        "Skipper baseline workspace")
    _require_exact_binding(workspace.sample_sum, binding.sample_sum,
        "Skipper sample-sum workspace")
    return nothing
end

@inline function _require_frame_readout_workspace_binding(
    workspace::MultiReadFrameReadoutWorkspace, binding::NamedTuple)
    _require_exact_binding(workspace.reference_average,
        binding.reference_average, "multi-read reference workspace")
    _require_exact_binding(workspace.signal_average,
        binding.signal_average, "multi-read signal workspace")
    _require_exact_binding(workspace.reference_cube,
        binding.reference_cube, "multi-read reference-cube workspace")
    _require_exact_binding(workspace.signal_cube,
        binding.signal_cube, "multi-read signal-cube workspace")
    return nothing
end

@inline function _require_frame_readout_workspace_binding(
    workspace::UpTheRampReadoutWorkspace, binding::NamedTuple)
    _require_exact_binding(workspace.slope, binding.slope,
        "ramp slope workspace")
    _require_exact_binding(workspace.intercept, binding.intercept,
        "ramp intercept workspace")
    _require_exact_binding(workspace.integrated, binding.integrated,
        "ramp integrated workspace")
    _require_exact_binding(workspace.cube, binding.cube,
        "ramp cube workspace")
    return nothing
end

@inline _require_frame_readout_product_binding(
    ::NoFrameReadoutProducts, ::NamedTuple{()}) = nothing

@inline function _require_frame_readout_product_binding(
    products::SkipperReadoutProducts, binding::NamedTuple)
    _require_exact_binding(products.mean_frame, binding.mean_frame,
        "Skipper mean-frame product")
    _require_exact_binding(products.sample_count, binding.sample_count,
        "Skipper sample-count product")
    return nothing
end

@inline function _require_frame_readout_product_binding(
    products::SampledFrameReadoutProducts, binding::NamedTuple)
    _require_exact_binding(products.reference_frame, binding.reference_frame,
        "sampled reference-frame product")
    _require_exact_binding(products.signal_frame, binding.signal_frame,
        "sampled signal-frame product")
    _require_exact_binding(products.read_cube, binding.read_cube,
        "sampled read-cube product")
    return nothing
end

@inline function _require_frame_readout_product_binding(
    products::MultiReadFrameReadoutProducts, binding::NamedTuple)
    _require_exact_binding(products.reference_frame, binding.reference_frame,
        "multi-read reference-frame product")
    _require_exact_binding(products.signal_frame, binding.signal_frame,
        "multi-read signal-frame product")
    _require_exact_binding(products.combined_frame, binding.combined_frame,
        "multi-read combined-frame product")
    _require_exact_binding(products.reference_cube, binding.reference_cube,
        "multi-read reference-cube product")
    _require_exact_binding(products.signal_cube, binding.signal_cube,
        "multi-read signal-cube product")
    _require_exact_binding(products.read_cube, binding.read_cube,
        "multi-read read-cube product")
    _require_exact_binding(products.read_offsets_s, binding.read_offsets_s,
        "multi-read duration product")
    return nothing
end

@inline function _require_frame_readout_product_binding(
    products::UpTheRampReadoutProducts, binding::NamedTuple)
    _require_exact_binding(products.slope_frame, binding.slope_frame,
        "ramp slope-frame product")
    _require_exact_binding(products.intercept_frame, binding.intercept_frame,
        "ramp intercept-frame product")
    _require_exact_binding(products.integrated_frame,
        binding.integrated_frame, "ramp integrated-frame product")
    _require_exact_binding(products.read_cube, binding.read_cube,
        "ramp read-cube product")
    _require_exact_binding(products.read_offsets_s, binding.read_offsets_s,
        "ramp time product")
    return nothing
end

@inline function _require_detector_state_binding(state::DetectorState,
    binding::NamedTuple)
    _require_exact_binding(state.accum_buffer, binding.accum_buffer,
        "detector accumulation state")
    _require_exact_binding(state.latent_buffer, binding.latent_buffer,
        "detector persistence state")
    _require_exact_binding(state.thermal_state, binding.thermal_state,
        "detector thermal state")
    return nothing
end

@inline function _require_detector_workspace_binding(
    workspace::DetectorWorkspace, binding::NamedTuple)
    _require_exact_binding(workspace.presampling_buffer,
        binding.presampling_buffer, "detector presampling workspace")
    _require_exact_binding(workspace.presampling_scratch,
        binding.presampling_scratch, "detector presampling scratch")
    _require_exact_binding(workspace.response_buffer, binding.response_buffer,
        "detector response workspace")
    _require_exact_binding(workspace.bin_buffer, binding.bin_buffer,
        "detector binning workspace")
    _require_exact_binding(workspace.temporal_buffer, binding.temporal_buffer,
        "detector temporal workspace")
    _require_exact_binding(workspace.noise_buffer, binding.noise_buffer,
        "detector noise workspace")
    _require_exact_binding(workspace.noise_buffer_host,
        binding.noise_buffer_host, "detector host-noise workspace")
    _require_exact_binding(workspace.batched_buffer_host,
        binding.batched_buffer_host, "detector batch workspace")
    _require_exact_binding(workspace.output_buffer_host,
        binding.output_buffer_host, "detector host-output workspace")
    _require_exact_binding(workspace.readout, binding.readout,
        "detector readout workspace owner")
    _require_frame_readout_workspace_binding(workspace.readout,
        binding.readout_binding)
    return nothing
end

@inline function _require_detector_product_binding(
    products::DetectorProducts, binding::NamedTuple)
    _require_exact_binding(products.frame, binding.frame,
        "detector frame product")
    _require_exact_binding(products.output_buffer, binding.output_buffer,
        "detector converted-output product")
    _require_exact_binding(products.readout, binding.readout,
        "detector readout-product owner")
    _require_frame_readout_product_binding(products.readout,
        binding.readout_binding)
    return nothing
end

@inline _require_no_alias(::Any, ::Any, ::AbstractString) = nothing

@inline function _require_no_alias(left::AbstractArray,
    right::AbstractArray, label::AbstractString)
    Base.mightalias(left, right) && throw(InvalidConfiguration(
        "$label must not alias"))
    return nothing
end

@inline _require_no_alias_against(::Any, ::Tuple{},
    ::AbstractString) = nothing

@inline function _require_no_alias_against(value, others::Tuple,
    label::AbstractString)
    _require_no_alias(value, first(others), label)
    return _require_no_alias_against(value, Base.tail(others), label)
end

@inline _require_no_internal_aliases(::Tuple{}, ::AbstractString) = nothing

@inline function _require_no_internal_aliases(members::Tuple,
    label::AbstractString)
    _require_no_alias_against(first(members), Base.tail(members), label)
    return _require_no_internal_aliases(Base.tail(members), label)
end

@inline _require_no_cross_aliases(::Tuple{}, ::Tuple,
    ::AbstractString) = nothing

@inline function _require_no_cross_aliases(left::Tuple, right::Tuple,
    label::AbstractString)
    _require_no_alias_against(first(left), right, label)
    return _require_no_cross_aliases(Base.tail(left), right, label)
end

function _require_disjoint_detector_owners(det::Detector,
    map::IntensityMap)
    state_members = _detector_state_members(det.state)
    workspace_members = _detector_workspace_members(det.workspace)
    product_members = _detector_product_members(det.products)
    _require_no_internal_aliases(state_members, "detector state storage")
    _require_no_internal_aliases(workspace_members,
        "detector workspace storage")
    _require_no_internal_aliases(product_members,
        "detector product storage")
    _require_no_cross_aliases(state_members, workspace_members,
        "detector state and workspace storage")
    _require_no_cross_aliases(state_members, product_members,
        "detector state and product storage")
    _require_no_cross_aliases(workspace_members, product_members,
        "detector workspace and product storage")
    _require_no_alias_against(map.values, state_members,
        "detector input and state storage")
    _require_no_alias_against(map.values, workspace_members,
        "detector input and workspace storage")
    _require_no_alias_against(map.values, product_members,
        "detector input and product storage")
    return nothing
end

function _detector_preparation_candidate(det::Detector)
    state = deepcopy(det.state)
    workspace = deepcopy(det.workspace)
    products = deepcopy(det.products)
    return Detector{
        typeof(det.noise),typeof(det.params),typeof(state),typeof(workspace),
        typeof(products),typeof(det.background_flux),
        typeof(det.background_map),typeof(backend(det)),
    }(det.noise, det.params, state, workspace, products,
        det.background_flux, det.background_map)
end

@inline function _require_detector_runtime_owner_types(det::Detector,
    candidate::Detector)
    typeof(candidate.state) === typeof(det.state) &&
        typeof(candidate.workspace) === typeof(det.workspace) &&
        typeof(candidate.products) === typeof(det.products) || throw(
        InvalidConfiguration(
            "prepared detector runtime owners do not fit the detector's concrete storage contract"))
    return nothing
end

function _detector_acquisition_preparation_contract(det::Detector,
    map::IntensityMap;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    require_whole_capture_idle(det)
    metadata = validate_plane_storage(map.metadata, map.values;
        label="detector intensity input")
    _require_detector_acquisition_plane(metadata.kind)
    _require_detector_incoherent(metadata.coherence)
    _require_finite_nonnegative_intensity(map.values)

    T = eltype(det.products.frame)
    metadata.numeric_type === T || throw(InvalidConfiguration(
        "detector and intensity map must use the same numeric type"))
    typeof(backend(det)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "detector and intensity map must use the same array backend"))
    compute_device(det.products.frame) == metadata.device ||
        throw(InvalidConfiguration(
            "detector and intensity map must occupy the same compute device"))

    normalization_scale = _normalization_rate_scale(metadata.normalization,
        normalized_to_photon_rate, T)
    spatial_scale = _spatial_rate_scale(metadata.spatial_measure, metadata, T)
    rate_scale = normalization_scale * spatial_scale
    isfinite(rate_scale) && rate_scale > zero(T) || throw(
        InvalidConfiguration(
            "prepared detector photon-rate scaling is not representable"))
    _require_prepared_response_sampling(det.params.response_model,
        det.params.psf_sampling)
    quantum_efficiency = _prepared_quantum_efficiency(det,
        metadata.spectral, T)
    frame_shape = detector_frame_shape(det, metadata.dimensions)
    output_shape = detector_output_shape(det, metadata.dimensions)
    _require_detector_configuration(det, frame_shape)
    _require_disjoint_detector_owners(det, map)
    return (; metadata, rate_scale, quantum_efficiency, frame_shape,
        output_shape)
end

function _prepare_detached_detector_acquisition(det::Detector,
    map::IntensityMap;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    contract = _detector_acquisition_preparation_contract(det, map;
        normalized_to_photon_rate=normalized_to_photon_rate)
    candidate = _detector_preparation_candidate(det)
    prepare_detector_buffers!(candidate, contract.metadata.dimensions)
    _require_prepared_detector_storage(candidate, contract.frame_shape,
        contract.output_shape)
    prepare_frame_readout_state!(candidate.params.sensor, candidate)
    _require_prepared_detector_storage(candidate, contract.frame_shape,
        contract.output_shape)
    _require_disjoint_detector_owners(candidate, map)
    _require_detector_runtime_owner_types(det, candidate)

    plan = DetectorAcquisitionPlan(_DETECTOR_ACQUISITION_PLAN_TOKEN,
        det.params, contract.metadata, contract.metadata.dimensions,
        contract.frame_shape, contract.output_shape, contract.rate_scale,
        contract.quantum_efficiency)
    return PreparedDetectorAcquisition(
        _PREPARED_DETECTOR_ACQUISITION_TOKEN, candidate, map, plan,
        candidate.state, candidate.workspace, candidate.products,
        _detector_state_binding(candidate.state),
        _detector_workspace_binding(candidate.workspace),
        _detector_product_binding(candidate.products))
end

function _rebind_prepared_detector_acquisition(det::Detector,
    candidate::PreparedDetectorAcquisition)
    candidate_detector = detector_acquisition_detector(candidate)
    det.params === candidate.plan.detector_params || throw(
        InvalidConfiguration(
            "detached detector acquisition does not match detector parameters"))
    _require_detector_runtime_owner_types(det, candidate_detector)
    return PreparedDetectorAcquisition(
        _PREPARED_DETECTOR_ACQUISITION_TOKEN, det, candidate.input,
        candidate.plan, candidate.state, candidate.workspace,
        candidate.products, candidate.state_binding,
        candidate.workspace_binding, candidate.product_binding)
end

@inline function _commit_prepared_detector_acquisition!(
    prepared::PreparedDetectorAcquisition)
    det = detector_acquisition_detector(prepared)
    det.state = prepared.state
    det.workspace = prepared.workspace
    det.products = prepared.products
    return prepared
end

"""
    prepare_detector_acquisition(detector, map;
        normalized_to_photon_rate=nothing)

Validate and prepare acquisition of an `IntensityMap`. Photon-rate maps cannot
be rescaled. A dimensionless calibration/test map must supply an explicit
`normalized_to_photon_rate` factor. Spatial-density samples are converted to a
cell rate using the declared plane sampling; cell-integrated samples are used
directly. Preparation also validates that every current sample is finite and
nonnegative.

The result is an exact `PreparedDetectorAcquisition`. Repeated `capture!` calls
trust producer writes to its bound intensity storage; producers must preserve
the finite, nonnegative value contract. Replacing any bound state, workspace,
product, or input storage requires preparation of a new owner.
"""
function prepare_detector_acquisition(det::Detector, map::IntensityMap;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    candidate = _prepare_detached_detector_acquisition(det, map;
        normalized_to_photon_rate=normalized_to_photon_rate)
    prepared = _rebind_prepared_detector_acquisition(det, candidate)
    _commit_prepared_detector_acquisition!(prepared)
    return prepared
end

@inline function _require_prepared_detector_binding(
    prepared::PreparedDetectorAcquisition)
    det = prepared.detector
    plan = prepared.plan
    det.params === plan.detector_params || throw(InvalidConfiguration(
        "detector parameters changed after acquisition preparation"))
    det.state === prepared.state || throw(InvalidConfiguration(
        "detector state owner changed after acquisition preparation"))
    det.workspace === prepared.workspace || throw(InvalidConfiguration(
        "detector workspace owner changed after acquisition preparation"))
    det.products === prepared.products || throw(InvalidConfiguration(
        "detector product owner changed after acquisition preparation"))
    prepared.input.metadata === plan.input_metadata || throw(
        InvalidConfiguration(
            "detector input metadata changed after acquisition preparation"))
    size(prepared.input.values) == plan.input_shape || throw(
        DimensionMismatchError(
            "detector input dimensions changed after acquisition preparation"))
    size(det.products.frame) == plan.frame_shape || throw(
        DimensionMismatchError(
            "detector frame dimensions changed after acquisition preparation"))
    size(output_frame(det)) == plan.output_shape || throw(
        DimensionMismatchError(
            "detector output dimensions changed after acquisition preparation"))
    typeof(backend(det)) === typeof(plan.input_metadata.backend) || throw(
        InvalidConfiguration(
            "detector backend changed after acquisition preparation"))
    compute_device(det.products.frame) == plan.input_metadata.device || throw(
        InvalidConfiguration(
            "detector device changed after acquisition preparation"))
    compute_device(prepared.input.values) == plan.input_metadata.device ||
        throw(InvalidConfiguration(
            "detector input device changed after acquisition preparation"))
    _require_detector_state_binding(det.state, prepared.state_binding)
    _require_detector_workspace_binding(det.workspace,
        prepared.workspace_binding)
    _require_detector_product_binding(det.products, prepared.product_binding)
    return nothing
end

@inline function _require_prepared_detector_binding(det::Detector,
    prepared::PreparedDetectorAcquisition)
    det === prepared.detector || throw(InvalidConfiguration(
        "detector does not match its prepared acquisition owner"))
    return _require_prepared_detector_binding(prepared)
end

@inline function _require_prepared_acquisition(
    prepared::PreparedDetectorAcquisition)
    _require_prepared_detector_binding(prepared)
    return nothing
end

@inline function _require_prepared_acquisition(
    prepared::PreparedDetectorAcquisition, map::IntensityMap)
    map.metadata === prepared.input.metadata &&
        map.values === prepared.input.values || throw(InvalidConfiguration(
        "intensity product does not match its prepared acquisition owner"))
    return _require_prepared_acquisition(prepared)
end

@inline function _require_prepared_whole_acquisition(
    prepared::PreparedDetectorAcquisition)
    _require_prepared_acquisition(prepared)
    require_whole_capture_idle(prepared.detector)
    return nothing
end

function capture!(prepared::PreparedDetectorAcquisition,
    rng::AbstractRNG)
    _require_prepared_whole_acquisition(prepared)
    plan = prepared.plan
    return capture_with_quantum_efficiency!(prepared.detector,
        prepared.input.values, plan.quantum_efficiency * plan.rate_scale, rng)
end

function capture!(prepared::PreparedDetectorAcquisition;
    rng::AbstractRNG=Random.default_rng(),
    integration_duration::Union{Nothing,Real}=nothing)
    integration_duration === nothing && return capture!(prepared, rng)
    _require_prepared_acquisition(prepared)
    plan = prepared.plan
    return capture_incremental!(prepared.detector, prepared.input.values, rng,
        integration_duration, plan.quantum_efficiency * plan.rate_scale)
end

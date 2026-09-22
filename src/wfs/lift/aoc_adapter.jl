"""
    LiFTForwardModel(forward)

Adapt the run-immutable physical plan of a CPU `PreparedLiFTForward` to
`AdaptiveOpticsCalibration.PhaseRetrieval`. The adapter does not retain the
prepared forward owner's OPD input, output, or replaceable scratch. Each
AdaptiveOpticsCalibration LiFT workspace receives independent AOS physical
forward scratch.
"""
struct LiFTForwardModel{T<:AbstractFloat,P<:LiFTForwardPlan{T}} <:
       PhaseRetrieval.AbstractLiFTForwardModel{T}
    plan::P
end

struct LiFTForwardModelWorkspace{
    T<:AbstractFloat,
    W<:LiFTForwardWorkspace,
    O<:AbstractMatrix{T},
}
    forward::W
    opd::O
end

@inline _require_lift_calibration_backend(::CPUBackend) = nothing
function _require_lift_calibration_backend(::AbstractArrayBackend)
    throw(UnsupportedAlgorithm(
        "AdaptiveOpticsCalibration LiFT integration currently supports CPU LiFT forward models only",
    ))
end

function LiFTForwardModel(forward::PreparedLiFTForward)
    _require_lift_forward_owner(forward)
    _require_lift_calibration_backend(forward.backend)
    plan = forward.plan
    T = eltype(plan.pupil_amplitude)
    return LiFTForwardModel{T,typeof(plan)}(plan)
end

@inline PhaseRetrieval.coefficient_count(model::LiFTForwardModel) =
    size(model.plan.basis, 3)

@inline PhaseRetrieval.observation_axes(model::LiFTForwardModel) =
    map(Base.OneTo, model.plan.observation_contract.rate_metadata.dimensions)

function PhaseRetrieval.allocate_model_workspace(model::LiFTForwardModel)
    plan = model.plan
    forward = _allocate_lift_forward_workspace(plan.pupil_amplitude,
        size(plan.pupil_amplitude, 1), plan.focal_resolution,
        plan.zero_padding, plan.object_kernel, plan.mapping)
    opd = similar(plan.pupil_amplitude)
    return LiFTForwardModelWorkspace(forward, opd)
end

function PhaseRetrieval.allocate_photon_rate(model::LiFTForwardModel{T}) where {T}
    plan = model.plan
    dimensions = plan.observation_contract.rate_metadata.dimensions
    return similar(plan.pupil_amplitude, T, dimensions...)
end

@inline function _aos_lift_model_workspace!(model::LiFTForwardModel,
    workspace::LiFTForwardModelWorkspace)
    _require_lift_forward_workspace(model.plan, workspace.forward)
    _require_lift_forward_input(model.plan, workspace.opd)
    (_lift_mightalias_any(workspace.opd,
        _lift_forward_plan_arrays(model.plan)) ||
        _lift_mightalias_forward_workspace(workspace.opd,
            workspace.forward)) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT OPD workspace must not alias physical-plan or forward scratch storage"))
    return workspace
end

@inline function _aos_lift_model_opd!(workspace::LiFTForwardModelWorkspace,
    model::LiFTForwardModel{T}, coefficients::AbstractVector{T}) where {T<:AbstractFloat}
    Base.require_one_based_indexing(coefficients)
    length(coefficients) == PhaseRetrieval.coefficient_count(model) || throw(
        DimensionMismatchError("AdaptiveOpticsCalibration LiFT coefficients must cover the complete AOS modal basis"))
    plan_arrays = _lift_forward_plan_arrays(model.plan)
    _lift_mightalias_any(coefficients, plan_arrays) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT coefficients must not alias physical-plan storage"))
    (Base.mightalias(coefficients, workspace.opd) ||
        _lift_mightalias_forward_workspace(coefficients,
            workspace.forward)) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT coefficients must not alias physical-model scratch"))
    plan = model.plan
    lift_basis_expansion!(workspace.opd, plan.basis, coefficients,
        plan.pupil_mask)
    @. workspace.opd += plan.diversity_opd
    return workspace.opd
end

@inline function _require_aos_lift_rate_output(model::LiFTForwardModel{T},
    workspace::LiFTForwardModelWorkspace,
    out::AbstractMatrix{T}, coefficients::AbstractVector{T}) where {T<:AbstractFloat}
    axes(out) == PhaseRetrieval.observation_axes(model) || throw(
        DimensionMismatchError("AdaptiveOpticsCalibration LiFT photon-rate output axes must match the AOS forward model"))
    eltype(out) === eltype(model.plan.pupil_amplitude) || throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT photon-rate output must use the AOS forward numeric type"))
    typeof(backend(out)) === typeof(backend(model.plan.pupil_amplitude)) || throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT photon-rate output must use the AOS forward backend"))
    compute_device(out) == compute_device(model.plan.pupil_amplitude) || throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT photon-rate output must occupy the AOS forward compute device"))
    (Base.mightalias(out, workspace.opd) ||
        _lift_mightalias_forward_workspace(out, workspace.forward)) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT photon-rate output must not alias model scratch"))
    _lift_mightalias_any(out, _lift_forward_plan_arrays(model.plan)) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT photon-rate output must not alias physical-plan storage"))
    Base.mightalias(out, coefficients) && throw(InvalidConfiguration(
        "AdaptiveOpticsCalibration LiFT photon-rate output must not alias coefficients"))
    return out
end

function PhaseRetrieval.predict_photon_rate!(out::AbstractMatrix{T},
    model::LiFTForwardModel{T}, workspace::LiFTForwardModelWorkspace,
    coefficients::AbstractVector{T}) where {T<:AbstractFloat}
    _aos_lift_model_workspace!(model, workspace)
    _require_aos_lift_rate_output(model, workspace, out, coefficients)
    opd = _aos_lift_model_opd!(workspace, model, coefficients)
    rate = _lift_rate_values_from_opd!(model.plan, workspace.forward, opd)
    copyto!(out, rate)
    return out
end

function PhaseRetrieval.analytic_photon_rate_jacobian!(jacobian::AbstractMatrix{T},
    model::LiFTForwardModel{T}, workspace::LiFTForwardModelWorkspace,
    coefficients::AbstractVector{T}, mode_indices::AbstractVector{Int}) where {T<:AbstractFloat}
    _aos_lift_model_workspace!(model, workspace)
    Base.require_one_based_indexing(jacobian, coefficients, mode_indices)
    expected_axes = (
        Base.OneTo(prod(model.plan.observation_contract.rate_metadata.dimensions)),
        Base.OneTo(length(mode_indices)),
    )
    axes(jacobian) == expected_axes || throw(DimensionMismatchError(
        "AdaptiveOpticsCalibration LiFT Jacobian dimensions must match the AOS observation and selected modes"))
    typeof(backend(jacobian)) === typeof(backend(model.plan.pupil_amplitude)) ||
        throw(InvalidConfiguration(
            "AdaptiveOpticsCalibration LiFT Jacobian must use the AOS forward backend"))
    compute_device(jacobian) == compute_device(model.plan.pupil_amplitude) ||
        throw(InvalidConfiguration(
            "AdaptiveOpticsCalibration LiFT Jacobian must occupy the AOS forward compute device"))
    plan_arrays = _lift_forward_plan_arrays(model.plan)
    _lift_mightalias_any(jacobian, plan_arrays) && throw(InvalidConfiguration(
        "AdaptiveOpticsCalibration LiFT Jacobian must not alias physical-plan storage"))
    (Base.mightalias(jacobian, workspace.opd) ||
        _lift_mightalias_forward_workspace(jacobian,
            workspace.forward)) && throw(
        InvalidConfiguration("AdaptiveOpticsCalibration LiFT Jacobian must not alias physical-model scratch"))
    Base.mightalias(jacobian, coefficients) && throw(InvalidConfiguration(
        "AdaptiveOpticsCalibration LiFT Jacobian must not alias coefficients"))
    coefficient_count = PhaseRetrieval.coefficient_count(model)
    @inbounds for mode_index in mode_indices
        1 <= mode_index <= coefficient_count || throw(DimensionMismatchError(
            "AdaptiveOpticsCalibration LiFT mode indices must select the AOS modal basis"))
    end
    _aos_lift_model_opd!(workspace, model, coefficients)
    plan = model.plan
    forward = workspace.forward
    amplitude_scale = sqrt(plan.photon_irradiance * plan.pupil_cell_area_m2)
    @. forward.amplitude_buffer = plan.pupil_amplitude * amplitude_scale
    oversampling = focal_field_from_opd!(forward.focal_buffer, plan, forward,
        forward.amplitude_buffer, workspace.opd)
    conjugate_field!(forward.conjugate_field_buffer, forward.focal_buffer)
    wavenumber = T(2 * pi) / plan.wavelength_m
    @inbounds for (column, mode_index) in enumerate(mode_indices)
        lift_scaled_basis_mode!(forward.amplitude_buffer, plan.pupil_amplitude,
            plan.basis, mode_index, amplitude_scale)
        focal_field_from_opd!(forward.mode_buffer, plan, forward,
            forward.amplitude_buffer, workspace.opd)
        field_derivative!(forward.optical_rate_buffer, forward.mode_buffer,
            forward.conjugate_field_buffer, oversampling, T(2) * wavenumber,
            forward.field_scratch)
        maybe_object_convolve!(plan, forward, forward.optical_rate_buffer)
        derivative = _apply_lift_mapping!(plan.mapping, forward,
            forward.optical_rate_buffer)
        lift_copy_column!(jacobian, column, derivative)
    end
    return jacobian
end

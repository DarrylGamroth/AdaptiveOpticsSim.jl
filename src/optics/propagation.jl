"""Nominal interface for a prepared propagation plan/workspace owner."""
abstract type AbstractPropagationModel end

"""
    AbstractPropagationPlan

Nominal interface for a run-immutable optical-propagation contract. Concrete
plans provide `propagation_input_metadata` and `propagation_output_metadata`.
They may own immutable backend-resident coefficients, but never FFT scratch or
a scratch-owning, stream-bound, or otherwise non-reentrant FFT handle.
"""
abstract type AbstractPropagationPlan end

struct FraunhoferPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_rad::T
end

"""Run-immutable Fraunhofer sampling and optical-plane contract."""
struct FraunhoferPropagationPlan{
    P<:FraunhoferPropagationParams,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
} <: AbstractPropagationPlan
    params::P
    input_metadata::I
    output_metadata::O
end

"""Replaceable single-writer scratch and FFT resource for Fraunhofer propagation."""
struct FraunhoferPropagationWorkspace{C<:AbstractMatrix,P}
    scratch::C
    fft_plan::P
end

"""Prepared Fraunhofer propagation plan/workspace pair."""
struct FraunhoferPropagation{
    P<:FraunhoferPropagationPlan,
    W<:FraunhoferPropagationWorkspace,
} <: AbstractPropagationModel
    plan::P
    workspace::W

    function FraunhoferPropagation(plan::P, workspace::W) where {
        P<:FraunhoferPropagationPlan,
        W<:FraunhoferPropagationWorkspace,
    }
        _require_propagation_plan_workspace(plan, workspace)
        return new{P,W}(plan, workspace)
    end
end

struct FresnelPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_m::T
    distance_m::T
end

"""Run-immutable Fresnel contract and backend-resident transfer operator."""
struct FresnelPropagationPlan{
    P<:FresnelPropagationParams,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
    M<:AbstractMatrix,
} <: AbstractPropagationPlan
    params::P
    input_metadata::I
    output_metadata::O
    transfer::M
end

"""Replaceable single-writer scratch and FFT resources for Fresnel propagation."""
struct FresnelPropagationWorkspace{C<:AbstractMatrix,Pf,Pi}
    spectrum::C
    propagated::C
    fft_plan::Pf
    ifft_plan::Pi
end

"""Prepared Fresnel propagation plan/workspace pair."""
struct FresnelPropagation{
    P<:FresnelPropagationPlan,
    W<:FresnelPropagationWorkspace,
} <: AbstractPropagationModel
    plan::P
    workspace::W

    function FresnelPropagation(plan::P, workspace::W) where {
        P<:FresnelPropagationPlan,
        W<:FresnelPropagationWorkspace,
    }
        _require_propagation_plan_workspace(plan, workspace)
        return new{P,W}(plan, workspace)
    end
end

"""Return the run-immutable plan held by a prepared propagation owner."""
@inline propagation_plan(model::AbstractPropagationModel) = model.plan

"""Return the single-writer workspace held by a prepared propagation owner."""
@inline propagation_workspace(model::AbstractPropagationModel) =
    model.workspace

"""Return the required input-plane metadata for a propagation plan or owner."""
@inline propagation_input_metadata(plan::AbstractPropagationPlan) =
    plan.input_metadata

"""Return the produced output-plane metadata for a propagation plan or owner."""
@inline propagation_output_metadata(plan::AbstractPropagationPlan) =
    plan.output_metadata
@inline propagation_input_metadata(model::AbstractPropagationModel) =
    propagation_input_metadata(propagation_plan(model))
@inline propagation_output_metadata(model::AbstractPropagationModel) =
    propagation_output_metadata(propagation_plan(model))

@kernel function complex_scale_copy_kernel!(out, input, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds out[i, j] = input[i, j] * scale
    end
end

@kernel function complex_hadamard_kernel!(out, weights, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds out[i, j] *= weights[i, j]
    end
end

@kernel function fresnel_transfer_kernel!(transfer, freqs, coeff, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds begin
            fx = freqs[i]
            fy = freqs[j]
            transfer[i, j] = cis(coeff * (fx * fx + fy * fy))
        end
    end
end

@inline function _fraunhofer_output_sampling(field::ElectricField)
    T = real(eltype(field.values))
    n = field.metadata.dimensions[1]
    input_sampling = field.metadata.sampling[1]
    return T(electric_field_wavelength(field) / (n * input_sampling))
end

function _fraunhofer_output_metadata(
    field::ElectricField,
    storage::AbstractMatrix,
    output_sampling,
    coherence::AbstractCombinationPolicy,
)
    return OpticalPlaneMetadata(
        FocalPlane(),
        storage;
        coordinate_domain=AngularCoordinates(),
        sampling=(output_sampling, output_sampling),
        orientation=field.metadata.orientation,
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=coherence,
    )
end

function _require_propagation_array_contract(
    storage::AbstractMatrix,
    metadata::OpticalPlaneMetadata,
    label::AbstractString,
)
    size(storage) == metadata.dimensions || throw(DimensionMismatchError(
        "$label dimensions do not match the propagation plan"))
    eltype(storage) === metadata.numeric_type || throw(InvalidConfiguration(
        "$label numeric type does not match the propagation plan"))
    typeof(backend(storage)) === typeof(metadata.backend) || throw(
        InvalidConfiguration(
            "$label backend does not match the propagation plan"))
    compute_device(storage) == metadata.device || throw(
        InvalidConfiguration(
            "$label device does not match the propagation plan"))
    return storage
end

function _require_propagation_plan_workspace(
    plan::FraunhoferPropagationPlan,
    workspace::FraunhoferPropagationWorkspace,
)
    _require_propagation_array_contract(
        workspace.scratch,
        plan.output_metadata,
        "Fraunhofer workspace scratch",
    )
    return workspace
end

function _require_propagation_plan_workspace(
    plan::FresnelPropagationPlan,
    workspace::FresnelPropagationWorkspace,
)
    _require_propagation_array_contract(
        plan.transfer,
        plan.output_metadata,
        "Fresnel transfer coefficients",
    )
    _require_propagation_array_contract(
        workspace.spectrum,
        plan.input_metadata,
        "Fresnel spectrum workspace",
    )
    _require_propagation_array_contract(
        workspace.propagated,
        plan.output_metadata,
        "Fresnel propagated workspace",
    )
    Base.mightalias(workspace.spectrum, workspace.propagated) && throw(
        InvalidConfiguration(
            "Fresnel spectrum and propagated workspaces must not alias"))
    Base.mightalias(workspace.spectrum, plan.transfer) && throw(
        InvalidConfiguration(
            "Fresnel spectrum workspace and transfer coefficients must not alias"))
    Base.mightalias(workspace.propagated, plan.transfer) && throw(
        InvalidConfiguration(
            "Fresnel propagated workspace and transfer coefficients must not alias"))
    return workspace
end

function FraunhoferPropagation(field::ElectricField)
    require_centered_plane_geometry(field.metadata;
        label="Fraunhofer input ElectricField")
    require_metric_coordinates(field.metadata;
        label="Fraunhofer input ElectricField")
    T = real(eltype(field.values))
    n = field.metadata.dimensions[1]
    n == field.metadata.dimensions[2] || throw(DimensionMismatchError(
        "Fraunhofer propagation requires a square ElectricField"))
    input_sampling = field.metadata.sampling[1]
    field.metadata.sampling[2] == input_sampling ||
        throw(InvalidConfiguration(
            "Fraunhofer propagation currently requires equal axis sampling"))
    scratch = similar(field.values)
    fft_plan = plan_fft_backend!(scratch)
    wavelength_m = electric_field_wavelength(field)
    output_sampling = _fraunhofer_output_sampling(field)
    params = FraunhoferPropagationParams{T}(
        n,
        wavelength_m,
        input_sampling,
        output_sampling,
    )
    output_metadata = _fraunhofer_output_metadata(
        field,
        scratch,
        output_sampling,
        field.metadata.coherence,
    )
    plan = FraunhoferPropagationPlan(params, field.metadata, output_metadata)
    workspace = FraunhoferPropagationWorkspace(scratch, fft_plan)
    return FraunhoferPropagation(plan, workspace)
end

function FresnelPropagation(field::ElectricField; distance_m::Real,
    output_kind::AbstractOpticalPlaneKind=IntermediatePlane())
    require_centered_plane_geometry(field.metadata;
        label="Fresnel input ElectricField")
    require_metric_coordinates(field.metadata;
        label="Fresnel input ElectricField")
    T = real(eltype(field.values))
    n = field.metadata.dimensions[1]
    n == field.metadata.dimensions[2] || throw(DimensionMismatchError(
        "Fresnel propagation requires a square ElectricField"))
    input_sampling = field.metadata.sampling[1]
    field.metadata.sampling[2] == input_sampling ||
        throw(InvalidConfiguration(
            "Fresnel propagation currently requires equal axis sampling"))
    spectrum = similar(field.values)
    propagated = similar(field.values)
    freqs = similar(field.values, T, n)
    transfer = similar(field.values)
    fft_plan = plan_fft_backend!(spectrum)
    ifft_plan = plan_ifft_backend!(propagated)
    params = FresnelPropagationParams{T}(
        n,
        electric_field_wavelength(field),
        input_sampling,
        input_sampling,
        T(distance_m),
    )
    output_metadata = OpticalPlaneMetadata(output_kind, propagated;
        coordinate_domain=MetricCoordinates(),
        sampling=field.metadata.sampling,
        orientation=field.metadata.orientation,
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=field.metadata.coherence)
    coeff = -T(pi) * params.wavelength * params.distance_m
    fftfreq!(freqs, params.padded_resolution; d=params.input_sampling_m)
    build_fresnel_transfer!(execution_style(transfer), transfer, freqs,
        coeff)
    plan = FresnelPropagationPlan(params, field.metadata, output_metadata,
        transfer)
    workspace = FresnelPropagationWorkspace(spectrum, propagated, fft_plan,
        ifft_plan)
    return FresnelPropagation(plan, workspace)
end

function FraunhoferPropagationWorkspace(
    plan::FraunhoferPropagationPlan,
    prototype::AbstractMatrix,
)
    _require_propagation_array_contract(
        prototype, plan.input_metadata, "Fraunhofer workspace prototype")
    scratch = similar(prototype)
    workspace = FraunhoferPropagationWorkspace(
        scratch, plan_fft_backend!(scratch))
    return _require_propagation_plan_workspace(plan, workspace)
end

function FresnelPropagationWorkspace(
    plan::FresnelPropagationPlan,
    prototype::AbstractMatrix,
)
    _require_propagation_array_contract(
        prototype, plan.input_metadata, "Fresnel workspace prototype")
    spectrum = similar(prototype)
    propagated = similar(prototype)
    workspace = FresnelPropagationWorkspace(
        spectrum,
        propagated,
        plan_fft_backend!(spectrum),
        plan_ifft_backend!(propagated),
    )
    return _require_propagation_plan_workspace(plan, workspace)
end

function _propagation_size_error(resolution::Int)
    throw(DimensionMismatchError("propagated field size must match propagation resolution $resolution"))
end

function _require_plan_match(field::ElectricField, plan::AbstractPropagationPlan)
    field.metadata == propagation_input_metadata(plan) || throw(InvalidConfiguration(
        "ElectricField metadata does not match the prepared propagation input"))
    return plan
end

@inline _require_model_match(field::ElectricField,
    model::AbstractPropagationModel) =
    (_require_plan_match(field, propagation_plan(model)); model)

function _complex_scale_copy!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, input::AbstractMatrix{Complex{T}}, scale::T) where {T<:AbstractFloat}
    n, m = size(out)
    @inbounds for j in 1:m, i in 1:n
        out[i, j] = input[i, j] * scale
    end
    return out
end

function _complex_scale_copy!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, input::AbstractMatrix{Complex{T}}, scale::T) where {T<:AbstractFloat}
    launch_kernel!(style, complex_scale_copy_kernel!, out, input, scale, size(out, 1); ndrange=size(out))
    return out
end

function _complex_hadamard!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, weights::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n, m = size(out)
    @inbounds for j in 1:m, i in 1:n
        out[i, j] *= weights[i, j]
    end
    return out
end

function _complex_hadamard!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, weights::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel!(style, complex_hadamard_kernel!, out, weights, size(out, 1); ndrange=size(out))
    return out
end

function build_fresnel_transfer!(::ScalarCPUStyle, transfer::AbstractMatrix{Complex{T}}, freqs::AbstractVector{T}, coeff::T) where {T<:AbstractFloat}
    n = length(freqs)
    @inbounds for j in 1:n, i in 1:n
        fx = freqs[i]
        fy = freqs[j]
        transfer[i, j] = cis(coeff * (fx * fx + fy * fy))
    end
    return transfer
end

function build_fresnel_transfer!(style::AcceleratorStyle, transfer::AbstractMatrix{Complex{T}}, freqs::AbstractVector{T}, coeff::T) where {T<:AbstractFloat}
    launch_kernel!(style, fresnel_transfer_kernel!, transfer, freqs, coeff, length(freqs); ndrange=size(transfer))
    return transfer
end

function propagate_fraunhofer_field!(out::AbstractMatrix{Complex{T}}, input::AbstractMatrix{Complex{T}}, scratch::AbstractMatrix{Complex{T}}, fft_plan) where {T<:AbstractFloat}
    size(out) == size(input) == size(scratch) || _propagation_size_error(size(input, 1))
    copyto!(scratch, input)
    execute_fft_plan!(scratch, fft_plan)
    _complex_scale_copy!(execution_style(out), out, scratch, inv(T(size(input, 1))))
    return out
end

function propagate_fresnel_field!(out::AbstractMatrix{Complex{T}}, input::AbstractMatrix{Complex{T}},
    spectrum::AbstractMatrix{Complex{T}}, propagated::AbstractMatrix{Complex{T}}, transfer::AbstractMatrix{Complex{T}},
    fft_plan, ifft_plan) where {T<:AbstractFloat}
    size(out) == size(input) == size(spectrum) == size(propagated) == size(transfer) ||
        _propagation_size_error(size(input, 1))
    copyto!(spectrum, input)
    execute_fft_plan!(spectrum, fft_plan)
    _complex_hadamard!(execution_style(spectrum), spectrum, transfer)
    copyto!(propagated, spectrum)
    execute_fft_plan!(propagated, ifft_plan)
    copyto!(out, propagated)
    return out
end

function propagate_field!(out::AbstractMatrix{Complex{T}},
    field::ElectricField, plan::FraunhoferPropagationPlan,
    workspace::FraunhoferPropagationWorkspace) where {T<:AbstractFloat}
    _require_plan_match(field, plan)
    size(out) == size(field.values) ||
        _propagation_size_error(plan.params.padded_resolution)
    require_same_backend(out, field)
    _require_propagation_plan_workspace(plan, workspace)
    return propagate_fraunhofer_field!(out, field.values,
        workspace.scratch, workspace.fft_plan)
end

@inline propagate_field!(out::AbstractMatrix{Complex{T}},
    field::ElectricField, model::FraunhoferPropagation) where {
    T<:AbstractFloat,
} = propagate_field!(out, field, model.plan, model.workspace)

function propagate_field!(out::ElectricField, field::ElectricField,
    model::FraunhoferPropagation)
    out.metadata == model.plan.output_metadata || throw(InvalidConfiguration(
        "destination ElectricField metadata does not match the prepared Fraunhofer output"))
    propagate_field!(out.values, field, model)
    return out
end

function propagate_field!(out::AbstractMatrix{Complex{T}},
    field::ElectricField, plan::FresnelPropagationPlan,
    workspace::FresnelPropagationWorkspace) where {T<:AbstractFloat}
    _require_plan_match(field, plan)
    size(out) == size(field.values) ||
        _propagation_size_error(plan.params.padded_resolution)
    require_same_backend(out, field)
    _require_propagation_plan_workspace(plan, workspace)
    return propagate_fresnel_field!(out, field.values,
        workspace.spectrum, workspace.propagated,
        plan.transfer, workspace.fft_plan, workspace.ifft_plan)
end

@inline propagate_field!(out::AbstractMatrix{Complex{T}},
    field::ElectricField, model::FresnelPropagation) where {
    T<:AbstractFloat,
} = propagate_field!(out, field, model.plan, model.workspace)

function propagate_field!(out::ElectricField, field::ElectricField,
    model::FresnelPropagation)
    out.metadata == model.plan.output_metadata || throw(InvalidConfiguration(
        "destination ElectricField metadata does not match the prepared Fresnel output"))
    propagate_field!(out.values, field, model)
    return out
end

function propagation_output(field::ElectricField,
    model::AbstractPropagationModel)
    _require_model_match(field, model)
    values = similar(field.values)
    return ElectricField(propagation_output_metadata(model), values)
end

function IntensityMap(field::ElectricField,
    model::FraunhoferPropagation)
    _require_model_match(field, model)
    _require_coherent_field(field.metadata.coherence)
    T = real(eltype(field.values))
    output_metadata = propagation_output_metadata(model)
    values = similar(field.values, T, output_metadata.dimensions...)
    fill!(values, zero(T))
    metadata = OpticalPlaneMetadata(output_metadata.kind, values;
        coordinate_domain=output_metadata.coordinate_domain,
        sampling=output_metadata.sampling,
        origin=output_metadata.origin,
        centering=output_metadata.centering,
        orientation=output_metadata.orientation,
        spectral=output_metadata.spectral,
        normalization=output_metadata.normalization,
        spatial_measure=output_metadata.spatial_measure,
        coherence=IncoherentIntensityAddition(),
        device=output_metadata.device)
    return IntensityMap(metadata, values)
end

function fraunhofer_intensity_from_field!(out::AbstractMatrix{T},
    field::ElectricField, plan::FraunhoferPropagationPlan,
    workspace::FraunhoferPropagationWorkspace) where {T<:AbstractFloat}
    _require_plan_match(field, plan)
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("Fraunhofer intensity output must match ElectricField size"))
    require_same_backend(out, field)
    _require_propagation_plan_workspace(plan, workspace)
    propagate_fraunhofer_field!(workspace.scratch, field.values,
        workspace.scratch, workspace.fft_plan)
    _intensity!(execution_style(out), out, workspace.scratch)
    return out
end

@inline fraunhofer_intensity_from_field!(out::AbstractMatrix{T},
    field::ElectricField, model::FraunhoferPropagation) where {
    T<:AbstractFloat,
} = fraunhofer_intensity_from_field!(out, field, model.plan,
    model.workspace)

function fraunhofer_intensity_from_field!(out::IntensityMap,
    field::ElectricField, model::FraunhoferPropagation)
    output_metadata = propagation_output_metadata(model)
    require_same_plane_grid(out.metadata, output_metadata;
        label="Fraunhofer intensity destination",
        require_numeric_type=false)
    require_compatible_radiometry(out.metadata, output_metadata;
        label="Fraunhofer intensity destination")
    _require_incoherent_policy(out.metadata.coherence,
        "Fraunhofer intensity destination")
    fraunhofer_intensity_from_field!(out.values, field, model)
    return out
end

function fraunhofer_intensity_stack!(intensity_stack::AbstractArray{T,3}, field_stack::AbstractArray{Complex{T},3}, fft_stack_plan) where {T<:AbstractFloat}
    execute_fft_plan!(field_stack, fft_stack_plan)
    intensity_scale = inv(T(size(field_stack, 1)) * T(size(field_stack, 2)))
    @. intensity_stack = abs2(field_stack) * intensity_scale
    return intensity_stack
end

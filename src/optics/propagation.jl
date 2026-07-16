abstract type AbstractPropagationModel end

struct FraunhoferPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_rad::T
end

struct FraunhoferPropagationState{C<:AbstractMatrix,P}
    scratch::C
    fft_plan::P
end

struct FraunhoferPropagation{
    P<:FraunhoferPropagationParams,
    S<:FraunhoferPropagationState,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
} <: AbstractPropagationModel
    params::P
    state::S
    input_metadata::I
    output_metadata::O
end

struct FresnelPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_m::T
    distance_m::T
end

struct FresnelPropagationState{C<:AbstractMatrix,V<:AbstractVector,M<:AbstractMatrix,Pf,Pi}
    spectrum::C
    propagated::C
    freqs::V
    transfer::M
    fft_plan::Pf
    ifft_plan::Pi
end

struct FresnelPropagation{
    P<:FresnelPropagationParams,
    S<:FresnelPropagationState,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
} <: AbstractPropagationModel
    params::P
    state::S
    input_metadata::I
    output_metadata::O
end

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
    output_sampling = T(wavelength_m / (n * input_sampling))
    params = FraunhoferPropagationParams{T}(
        n,
        wavelength_m,
        input_sampling,
        output_sampling,
    )
    state = FraunhoferPropagationState{typeof(scratch), typeof(fft_plan)}(scratch, fft_plan)
    output_metadata = OpticalPlaneMetadata(FocalPlane(), scratch;
        coordinate_domain=AngularCoordinates(),
        sampling=(output_sampling, output_sampling),
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=field.metadata.coherence)
    return FraunhoferPropagation(params, state, field.metadata,
        output_metadata)
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
    state = FresnelPropagationState{typeof(spectrum), typeof(freqs), typeof(transfer), typeof(fft_plan), typeof(ifft_plan)}(
        spectrum,
        propagated,
        freqs,
        transfer,
        fft_plan,
        ifft_plan,
    )
    output_metadata = OpticalPlaneMetadata(output_kind, propagated;
        coordinate_domain=MetricCoordinates(),
        sampling=field.metadata.sampling,
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=field.metadata.coherence)
    model = FresnelPropagation(params, state, field.metadata,
        output_metadata)
    build_fresnel_transfer!(model)
    return model
end

function _propagation_size_error(resolution::Int)
    throw(DimensionMismatchError("propagated field size must match propagation resolution $resolution"))
end

function _require_model_match(field::ElectricField, model::AbstractPropagationModel)
    field.metadata == model.input_metadata || throw(InvalidConfiguration(
        "ElectricField metadata does not match the prepared propagation input"))
    return model
end

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

function build_fresnel_transfer!(model::FresnelPropagation)
    fftfreq!(model.state.freqs, model.params.padded_resolution; d=model.params.input_sampling_m)
    coeff = -eltype(model.state.freqs)(pi) * model.params.wavelength * model.params.distance_m
    return build_fresnel_transfer!(execution_style(model.state.transfer), model.state.transfer, model.state.freqs, coeff)
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
    field::ElectricField, model::FraunhoferPropagation) where {T<:AbstractFloat}
    _require_model_match(field, model)
    size(out) == size(field.values) ||
        _propagation_size_error(model.params.padded_resolution)
    require_same_backend(out, field)
    return propagate_fraunhofer_field!(out, field.values,
        model.state.scratch, model.state.fft_plan)
end

function propagate_field!(out::ElectricField, field::ElectricField,
    model::FraunhoferPropagation)
    out.metadata == model.output_metadata || throw(InvalidConfiguration(
        "destination ElectricField metadata does not match the prepared Fraunhofer output"))
    propagate_field!(out.values, field, model)
    return out
end

function propagate_field!(out::AbstractMatrix{Complex{T}},
    field::ElectricField, model::FresnelPropagation) where {T<:AbstractFloat}
    _require_model_match(field, model)
    size(out) == size(field.values) ||
        _propagation_size_error(model.params.padded_resolution)
    require_same_backend(out, field)
    return propagate_fresnel_field!(out, field.values,
        model.state.spectrum, model.state.propagated,
        model.state.transfer, model.state.fft_plan, model.state.ifft_plan)
end

function propagate_field!(out::ElectricField, field::ElectricField,
    model::FresnelPropagation)
    out.metadata == model.output_metadata || throw(InvalidConfiguration(
        "destination ElectricField metadata does not match the prepared Fresnel output"))
    propagate_field!(out.values, field, model)
    return out
end

function propagation_output(field::ElectricField,
    model::AbstractPropagationModel)
    _require_model_match(field, model)
    values = similar(field.values)
    return ElectricField(model.output_metadata, values)
end

function IntensityMap(field::ElectricField,
    model::FraunhoferPropagation)
    _require_model_match(field, model)
    T = real(eltype(field.values))
    values = similar(field.values, T, model.output_metadata.dimensions...)
    metadata = OpticalPlaneMetadata(model.output_metadata.kind, values;
        coordinate_domain=model.output_metadata.coordinate_domain,
        sampling=model.output_metadata.sampling,
        origin=model.output_metadata.origin,
        centering=model.output_metadata.centering,
        orientation=model.output_metadata.orientation,
        spectral=model.output_metadata.spectral,
        normalization=model.output_metadata.normalization,
        spatial_measure=model.output_metadata.spatial_measure,
        coherence=model.output_metadata.coherence,
        device=model.output_metadata.device)
    return IntensityMap(metadata, values)
end

function fraunhofer_intensity_from_field!(out::AbstractMatrix{T},
    field::ElectricField, model::FraunhoferPropagation) where {T<:AbstractFloat}
    _require_model_match(field, model)
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("Fraunhofer intensity output must match ElectricField size"))
    require_same_backend(out, field)
    propagate_fraunhofer_field!(model.state.scratch, field.values,
        model.state.scratch, model.state.fft_plan)
    _intensity!(execution_style(out), out, model.state.scratch)
    return out
end

function fraunhofer_intensity_from_field!(out::IntensityMap,
    field::ElectricField, model::FraunhoferPropagation)
    require_same_plane_grid(out.metadata, model.output_metadata;
        label="Fraunhofer intensity destination",
        require_numeric_type=false)
    fraunhofer_intensity_from_field!(out.values, field, model)
    return out
end

function fraunhofer_intensity_stack!(intensity_stack::AbstractArray{T,3}, field_stack::AbstractArray{Complex{T},3}, fft_stack_plan) where {T<:AbstractFloat}
    execute_fft_plan!(field_stack, fft_stack_plan)
    @. intensity_stack = abs2(field_stack)
    return intensity_stack
end

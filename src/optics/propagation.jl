abstract type AbstractPropagationModel end

struct FraunhoferPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_rad::T
end

mutable struct FraunhoferPropagationState{C<:AbstractMatrix,P}
    scratch::C
    fft_plan::P
end

struct FraunhoferPropagation{P<:FraunhoferPropagationParams,S<:FraunhoferPropagationState} <: AbstractPropagationModel
    params::P
    state::S
end

struct FresnelPropagationParams{T<:AbstractFloat}
    padded_resolution::Int
    wavelength::T
    input_sampling_m::T
    output_sampling_m::T
    distance_m::T
end

mutable struct FresnelPropagationState{C<:AbstractMatrix,V<:AbstractVector,M<:AbstractMatrix,Pf,Pi}
    spectrum::C
    propagated::C
    freqs::V
    transfer::M
    fft_plan::Pf
    ifft_plan::Pi
end

struct FresnelPropagation{P<:FresnelPropagationParams,S<:FresnelPropagationState} <: AbstractPropagationModel
    params::P
    state::S
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
    T = eltype(field.state.intensity)
    n = field.params.padded_resolution
    scratch = similar(field.state.field)
    fft_plan = plan_fft_backend!(scratch)
    params = FraunhoferPropagationParams{T}(
        n,
        field.params.wavelength,
        field.params.sampling_m,
        T(field.params.wavelength / (n * field.params.sampling_m)),
    )
    state = FraunhoferPropagationState{typeof(scratch), typeof(fft_plan)}(scratch, fft_plan)
    return FraunhoferPropagation(params, state)
end

function FresnelPropagation(field::ElectricField; distance_m::Real)
    T = eltype(field.state.intensity)
    n = field.params.padded_resolution
    spectrum = similar(field.state.field)
    propagated = similar(field.state.field)
    freqs = similar(field.state.intensity, n)
    transfer = similar(field.state.field)
    fft_plan = plan_fft_backend!(spectrum)
    ifft_plan = plan_ifft_backend!(propagated)
    params = FresnelPropagationParams{T}(
        n,
        field.params.wavelength,
        field.params.sampling_m,
        field.params.sampling_m,
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
    model = FresnelPropagation(params, state)
    build_fresnel_transfer!(model)
    return model
end

function _propagation_size_error(resolution::Int)
    throw(DimensionMismatchError("propagated field size must match propagation resolution $resolution"))
end

function _require_model_match(field::ElectricField, model::AbstractPropagationModel)
    field.params.padded_resolution == model.params.padded_resolution ||
        throw(DimensionMismatchError("ElectricField padded resolution must match propagation model resolution"))
    field.params.wavelength == model.params.wavelength ||
        throw(InvalidConfiguration("ElectricField wavelength must match propagation model wavelength"))
    field.params.sampling_m == model.params.input_sampling_m ||
        throw(InvalidConfiguration("ElectricField sampling must match propagation model sampling"))
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

function propagate_field!(out::AbstractMatrix{Complex{T}}, field::ElectricField, model::FraunhoferPropagation) where {T<:AbstractFloat}
    _require_model_match(field, model)
    size(out) == size(field.state.field) || _propagation_size_error(model.params.padded_resolution)
    return propagate_fraunhofer_field!(out, field.state.field, model.state.scratch, model.state.fft_plan)
end

function propagate_field!(field::ElectricField, model::FraunhoferPropagation)
    propagate_field!(field.state.field, field, model)
    return field
end

function propagate_field!(out::AbstractMatrix{Complex{T}}, field::ElectricField, model::FresnelPropagation) where {T<:AbstractFloat}
    _require_model_match(field, model)
    size(out) == size(field.state.field) || _propagation_size_error(model.params.padded_resolution)
    return propagate_fresnel_field!(out, field.state.field, model.state.spectrum, model.state.propagated,
        model.state.transfer, model.state.fft_plan, model.state.ifft_plan)
end

function propagate_field!(field::ElectricField, model::FresnelPropagation)
    propagate_field!(field.state.field, field, model)
    return field
end

function fraunhofer_intensity_from_field!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.state.field) ||
        throw(DimensionMismatchError("Fraunhofer intensity output must match ElectricField size"))
    propagate_fraunhofer_field!(field.state.fft_buffer, field.state.field, field.state.fft_buffer, field.state.fft_plan)
    @. out = abs2(field.state.fft_buffer)
    return out
end

function fraunhofer_intensity_stack!(intensity_stack::AbstractArray{T,3}, field_stack::AbstractArray{Complex{T},3}, fft_stack_plan) where {T<:AbstractFloat}
    execute_fft_plan!(field_stack, fft_stack_plan)
    @. intensity_stack = abs2(field_stack)
    return intensity_stack
end

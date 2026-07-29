"""
Immutable model and numerical-sampling parameters for a square regular
microlens array.

`n_lenslets` is the number of lenslets along each pupil axis, and
`pixel_scale_arcsec` is the requested detector-plane angular sampling in
arcseconds per pixel. A composed Shack–Hartmann front end supplies the physical
subaperture pitch. This model does not currently describe focal length, fill
factor, per-lenslet prescriptions, or manufacturing errors.
"""
struct MicrolensArrayParams{T<:AbstractFloat}
    n_lenslets::Int
    half_pixel_shift::Bool
    diffraction_padding::Int
    pixel_scale_arcsec::Union{T,Nothing}
    n_pix_subap::Union{Int,Nothing}
    shannon_sampling::Bool
end

"""
A regular angular-coordinate microlens-array model, independent of detector
acquisition and WFS estimation.
"""
struct MicrolensArray{P<:MicrolensArrayParams}
    params::P
end

"""Return the microlens-array component composed into an optical front end."""
function microlens_array end

function MicrolensArray(; n_lenslets::Integer,
    half_pixel_shift::Bool=false,
    diffraction_padding::Integer=2, pixel_scale_arcsec=nothing,
    n_pix_subap::Union{Nothing,Integer}=nothing,
    shannon_sampling::Bool=true,
    T::Type{<:AbstractFloat}=pixel_scale_arcsec === nothing ? Float64 :
        typeof(float(pixel_scale_arcsec)))
    lenslets = Int(n_lenslets)
    padding = Int(diffraction_padding)
    pixels = n_pix_subap === nothing ? nothing : Int(n_pix_subap)
    lenslets > 0 || throw(InvalidConfiguration(
        "n_lenslets must be positive"))
    padding > 0 || throw(InvalidConfiguration(
        "diffraction_padding must be positive"))
    if pixels !== nothing
        pixels > 0 || throw(InvalidConfiguration(
            "n_pix_subap must be positive"))
        iseven(pixels) || throw(InvalidConfiguration(
            "n_pix_subap must be even"))
    end
    scale_arcsec = pixel_scale_arcsec === nothing ? nothing :
        T(pixel_scale_arcsec)
    if scale_arcsec !== nothing
        isfinite(scale_arcsec) && scale_arcsec > zero(T) ||
            throw(InvalidConfiguration(
                "pixel_scale_arcsec must be finite and positive"))
    end
    params = MicrolensArrayParams{T}(lenslets, half_pixel_shift,
        padding, scale_arcsec, pixels, shannon_sampling)
    return MicrolensArray{typeof(params)}(params)
end

@inline microlens_numeric_type(
    ::MicrolensArray{<:MicrolensArrayParams{T}}) where {T} = T

"""
Backend-bound plans and preallocated scratch for microlens propagation.

This is prepared execution state, not part of the physical `MicrolensArray`
description. One optical front end owns one instance under a single-writer
contract.
"""
mutable struct PreparedMicrolensPropagation{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    CC<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RC3<:AbstractArray{T,3},
    RT3<:AbstractArray{T,3},
    RB<:AbstractMatrix{T},
    RS<:AbstractMatrix{T},
    RCS<:AbstractArray{T,3},
    RA<:AbstractArray{T,3},
    P,
    PS,
    Pi,
    PSi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3},
    V<:AbstractVector{T}}
    field::C
    phasor::C
    fft_buffer::C
    fft_stack::CC
    intensity::R
    intensity_stack::RC3
    intensity_tmp_stack::RT3
    temp::R
    bin_buffer::RB
    spot::RS
    sampled_spot_cube::RCS
    spot_cube_accum::RA
    fft_plan::P
    fft_stack_plan::PS
    ifft_plan::Pi
    ifft_stack_plan::PSi
    elongation_kernel::K
    lgs_kernel_fft::KF
    lgs_kernel_tag::UInt
    effective_padding::Int
    binning_pixel_scale::Int
    sampled_n_pix_subap::Int
    phasor_ratio::T
    fft_asterism_stack::CC
    fft_asterism_plan::PS
    asterism_capacity::Int
    amp_scales::V
    amp_scales_host::Vector{T}
    opd_to_cycles::V
    opd_to_cycles_host::Vector{T}
end

"""
    prepare_microlens_propagation(microlens_array, pupil_resolution;
        T=Float64, backend=CPUBackend())

Prepare backend-specific FFT plans and reusable scratch storage for an
independent `MicrolensArray`. This cold-path constructor performs no detector
acquisition or slope-estimator setup.
"""
function prepare_microlens_propagation(mla::MicrolensArray,
    pupil_resolution::Integer;
    T::Type{<:AbstractFloat}=microlens_numeric_type(mla),
    backend::AbstractArrayBackend=CPUBackend())
    resolution = Int(pupil_resolution)
    resolution > 0 || throw(InvalidConfiguration(
        "pupil_resolution must be positive"))
    resolution % mla.params.n_lenslets == 0 ||
        throw(InvalidConfiguration(
            "pupil_resolution must be divisible by n_lenslets"))
    selector = _resolve_backend_selector(backend)
    array_type = _resolve_array_backend(selector)
    return _prepare_microlens_propagation(array_type, T, mla,
        div(resolution, mla.params.n_lenslets))
end

function _prepare_microlens_propagation(backend, ::Type{T},
    mla::MicrolensArray, sub::Int) where {T<:AbstractFloat}
    n_lenslets = mla.params.n_lenslets
    pad = max(sub, sub * mla.params.diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    phasor = similar(field)
    fft_buffer = similar(field)
    fft_stack = backend{Complex{T}}(undef, pad, pad,
        n_lenslets * n_lenslets)
    intensity = backend{T}(undef, pad, pad)
    intensity_stack = backend{T}(undef, pad, pad,
        n_lenslets * n_lenslets)
    intensity_tmp_stack = similar(intensity_stack)
    temp = similar(intensity)
    bin_buffer = backend{T}(undef, sub, sub)
    spot = similar(bin_buffer)
    sampled_spot_cube = backend{T}(undef, n_lenslets * n_lenslets,
        sub, sub)
    spot_cube_accum = similar(sampled_spot_cube)
    fft_plan = plan_fft_backend!(fft_buffer)
    fft_stack_plan = plan_fft_backend!(fft_stack, (1, 2))
    ifft_plan = plan_ifft_backend!(fft_buffer)
    ifft_stack_plan = plan_ifft_backend!(fft_stack, (1, 2))
    fft_asterism_stack = similar(fft_stack)
    fft_asterism_plan = plan_fft_backend!(fft_asterism_stack, (1, 2))
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0, 0)
    amp_scales = backend{T}(undef, 1)
    amp_scales_host = Vector{T}(undef, 1)
    opd_to_cycles = backend{T}(undef, 1)
    opd_to_cycles_host = Vector{T}(undef, 1)
    propagation = PreparedMicrolensPropagation(field, phasor, fft_buffer,
        fft_stack, intensity, intensity_stack, intensity_tmp_stack, temp,
        bin_buffer, spot, sampled_spot_cube, spot_cube_accum, fft_plan,
        fft_stack_plan, ifft_plan, ifft_stack_plan, elongation_kernel,
        lgs_kernel_fft, UInt(0), mla.params.diffraction_padding, 1, sub,
        T(NaN), fft_asterism_stack, fft_asterism_plan, 1, amp_scales,
        amp_scales_host, opd_to_cycles, opd_to_cycles_host)
    return propagation
end

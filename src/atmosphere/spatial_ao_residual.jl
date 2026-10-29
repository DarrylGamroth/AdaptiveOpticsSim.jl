import ..Backends: enqueue_fft_plan!, plan_fft_backend!

"""
    SpatialAORejection(; control_cutoff_cycles_per_m, residual_amplitude)

Explicit idealized spatial AO rejection on each full periodic moving screen.
The real, even Fourier-amplitude transfer is zero at DC, `residual_amplitude`
for nonzero radial frequencies at or below `control_cutoff_cycles_per_m`, and
one above it. The PSD fraction below cutoff is `residual_amplitude^2`.

The cutoff must be positive and no greater than the screen's axis Nyquist
frequency (checked during preparation). The amplitude must lie in [0, 1].
This instantaneous model has no WFS noise, temporal servo, or facility fit.
"""
struct SpatialAORejection{T<:AbstractFloat}
    control_cutoff_cycles_per_m::T
    residual_amplitude::T

    function SpatialAORejection(; control_cutoff_cycles_per_m::Real,
            residual_amplitude::Real)
        control_cutoff_cycles_per_m isa Bool && throw(InvalidConfiguration(
            "spatial AO control cutoff must be a real frequency, not Bool"))
        residual_amplitude isa Bool && throw(InvalidConfiguration(
            "spatial AO residual amplitude must be a real fraction, not Bool"))
        cutoff, amplitude = promote(float(control_cutoff_cycles_per_m),
            float(residual_amplitude))
        isfinite(cutoff) && cutoff > zero(cutoff) || throw(InvalidConfiguration(
            "spatial AO control cutoff must be finite and positive"))
        isfinite(amplitude) && zero(amplitude) <= amplitude <= one(amplitude) ||
            throw(InvalidConfiguration("spatial AO residual amplitude must lie in [0, 1]"))
        return new{typeof(cutoff)}(cutoff, amplitude)
    end
end

struct _SpatialAORejectionWorkspace{A,F,I}
    buffer::A
    fft_plan::F
    ifft_plan::I
end

# Exact screen bindings and owned transfer survive reset; only the FFT buffer
# is replaceable scratch. Filtering never advances the atmosphere or its RNG.
struct _PreparedSpatialAORejection{R,A,S,W}
    rejection::R
    transfer::A
    screens::S
    workspace::W
end

function _spatial_ao_transfer(rejection::SpatialAORejection, n::Int,
        sampling_m::T) where {T<:AbstractFloat}
    rejection.control_cutoff_cycles_per_m <= inv(T(2) * sampling_m) ||
        throw(InvalidConfiguration("spatial AO control cutoff exceeds axis Nyquist"))
    cutoff = T(rejection.control_cutoff_cycles_per_m)
    amplitude = T(rejection.residual_amplitude)
    isfinite(cutoff) && zero(T) < cutoff <= inv(T(2) * sampling_m) ||
        throw(InvalidConfiguration(
            "spatial AO control cutoff must be representable and not exceed axis Nyquist"))
    isfinite(amplitude) && zero(T) <= amplitude <= one(T) ||
        throw(InvalidConfiguration("spatial AO amplitude must remain in [0, 1] at screen precision"))
    freqs = Vector{T}(undef, n)
    fftfreq!(freqs, n; d=sampling_m)
    transfer = Matrix{T}(undef, n, n)
    for j in axes(transfer, 2), i in axes(transfer, 1)
        transfer[i, j] = i == 1 && j == 1 ? zero(T) :
            hypot(freqs[i], freqs[j]) <= cutoff ? amplitude : one(T)
    end
    return transfer
end

_prepare_spatial_ao_rejection(::Nothing, ::MultiLayerAtmosphere) = nothing

function _prepare_spatial_ao_rejection(rejection::SpatialAORejection,
        atmosphere::MultiLayerAtmosphere)
    layers = atmosphere.layers
    screens = map(layer -> layer.generator.state.phase_rad, Tuple(layers))
    first_screen = first(screens)
    T = eltype(first_screen)
    n = size(first_screen, 1)
    sampling_m = first(layers).params.sampling_m
    device = compute_device(first_screen)
    for (layer, screen) in zip(layers, screens)
        size(screen) == (n, n) && eltype(screen) === T || throw(InvalidConfiguration(
            "spatial AO screens must have the same square grid and precision"))
        layer.params.sampling_m == sampling_m || throw(InvalidConfiguration(
            "spatial AO screens must have the same physical sampling"))
        compute_device(screen) == device || throw(InvalidConfiguration(
            "spatial AO screens must occupy the same compute device"))
    end
    host_transfer = _spatial_ao_transfer(rejection, n, sampling_m)
    transfer, workspace = _with_compute_device(device) do
        transfer = allocate_device_array(device, T, n, n)
        copyto!(transfer, host_transfer)
        buffer = allocate_device_array(device, Complex{T}, n, n)
        fft_plan = plan_fft_backend!(buffer)
        ifft_plan = plan_ifft_backend!(buffer)
        return transfer, _SpatialAORejectionWorkspace(buffer, fft_plan, ifft_plan)
    end
    return _PreparedSpatialAORejection(rejection, transfer, screens, workspace)
end

_apply_spatial_ao_rejection!(::Nothing) = nothing

function _apply_spatial_ao_rejection!(owner::_PreparedSpatialAORejection)
    workspace = owner.workspace
    buffer = workspace.buffer
    for screen in owner.screens
        buffer .= screen
        enqueue_fft_plan!(buffer, workspace.fft_plan)
        buffer .*= owner.transfer
        enqueue_fft_plan!(buffer, workspace.ifft_plan)
        screen .= real.(buffer)
    end
    return nothing
end

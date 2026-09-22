#
# Zernike wavefront sensing
#
# The maintained MVP follows the standard phase-shifting Zernike-wavefront-
# sensor optical chain:
#
# 1. embed the pupil field on a padded FFT grid
# 2. propagate to the focal plane
# 3. apply a circular phase-shifting spot of radius `spot_radius_lambda_over_d`
# 4. propagate back to the pupil plane
# 5. sample the re-imaged pupil intensity on a compact detector grid
# 6. hand the photon-arrival-rate frame to a detector or external estimator
#

@kernel function zernike_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end


"""Run-immutable numerical contract for phase-spot pupil-relay propagation."""
struct ZernikePropagationPlan{M<:ZernikePhaseSpot,T<:AbstractFloat}
    phase_spot::M
    pupil_resolution::Int
    pupil_samples::Int
    numeric_type::Type{T}
end

"""
Backend-bound FFT handles and replaceable single-writer scratch for phase-spot
pupil-relay propagation. No field is a caller-visible optical product.
"""
mutable struct ZernikePropagationWorkspace{
    T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
}
    field::C
    focal_field::C
    pupil_field::C
    phasor::C
    phase_mask::C
    pupil_intensity::R
    nominal_frame::R
    fft_plan::Pf
    ifft_plan::Pi
end

"""Exact plan/workspace owner for one Zernike propagation execution."""
struct PreparedZernikePropagation{
    P<:ZernikePropagationPlan,W<:ZernikePropagationWorkspace}
    plan::P
    workspace::W
end

@inline zernike_propagation_plan(
    propagation::PreparedZernikePropagation) = propagation.plan
@inline zernike_propagation_workspace(
    propagation::PreparedZernikePropagation) = propagation.workspace

"""Zernike phase spot and prepared re-imaged-pupil optical front end."""
struct ZernikeOpticalFrontEnd{M,P,S}
    phase_spot::M
    propagation::P
    binning::Int
    source::S
end

"""Run-immutable binning contract for the internal detector-facing frame."""
struct ZernikeAcquisitionPlan
    binning::Int
end

"""Caller-visible detector-facing photon-arrival-rate frame product."""
mutable struct ZernikeAcquisitionProducts{T<:AbstractFloat,
    R<:AbstractMatrix{T}}
    frame::R
end

"""
Internal convenience owner for frame binning and its rate product. Detector
response, exposure integration, and readout belong to the generic acquisition
stage.
"""
struct ZernikeDetectorAcquisition{P,PR}
    plan::P
    products::PR
end

"""
    ZernikeWFS

Diffractive Zernike wavefront sensor with a circular focal-plane phase spot.
It produces physical photon-arrival-rate frames; estimation belongs to the RTC.
"""
struct ZernikeWFS{F,A,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
end

@inline backend(::ZernikeWFS{F,A,B}) where {F,A,B} = B()

@inline zernike_acquisition_plan(wfs::ZernikeWFS) = wfs.acquisition.plan
@inline zernike_acquisition_products(wfs::ZernikeWFS) =
    wfs.acquisition.products
@inline zernike_propagation(wfs::ZernikeWFS) = wfs.front_end.propagation
@inline zernike_propagation_plan(wfs::ZernikeWFS) =
    zernike_propagation_plan(zernike_propagation(wfs))
@inline zernike_propagation_workspace(wfs::ZernikeWFS) =
    zernike_propagation_workspace(zernike_propagation(wfs))
@inline zernike_propagation_workspace(front_end::ZernikeOpticalFrontEnd) =
    zernike_propagation_workspace(front_end.propagation)

"""
    ZernikeWFS(tel; ...)

Construct a Zernike WFS using a focal-plane circular phase-shifting spot.

`pupil_samples` defines the nominal sampled pupil grid before optional `binning`
coarsens the final exported camera/signal frame.
"""
function ZernikeWFS(tel::Telescope; pupil_samples::Int,
    phase_shift_pi::Real=0.5,
    spot_radius_lambda_over_d::Real=1.0,
    diffraction_padding::Int=2,
    binning::Int=1,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    if tel.params.resolution % pupil_samples != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by pupil_samples"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if pupil_samples % binning != 0
        throw(InvalidConfiguration("pupil_samples must be divisible by binning"))
    end
    if diffraction_padding < 1
        throw(InvalidConfiguration("diffraction_padding must be >= 1"))
    end
    n_signal = div(pupil_samples, binning)
    pad = tel.params.resolution * diffraction_padding
    spot = ZernikePhaseSpot(T(phase_shift_pi),
        T(spot_radius_lambda_over_d), diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    phasor = similar(field)
    phase_mask = similar(field)
    pupil_intensity = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    nominal_frame = backend{T}(undef, pupil_samples, pupil_samples)
    camera_frame = backend{T}(undef, n_signal, n_signal)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    propagation_workspace = ZernikePropagationWorkspace{
        T,
        typeof(field),
        typeof(pupil_intensity),
        typeof(fft_plan),
        typeof(ifft_plan),
    }(
        field,
        focal_field,
        pupil_field,
        phasor,
        phase_mask,
        pupil_intensity,
        nominal_frame,
        fft_plan,
        ifft_plan,
    )
    propagation_plan = ZernikePropagationPlan(spot,
        tel.params.resolution, pupil_samples, T)
    propagation = PreparedZernikePropagation(propagation_plan,
        propagation_workspace)
    acquisition = ZernikeDetectorAcquisition(
        ZernikeAcquisitionPlan(binning),
        ZernikeAcquisitionProducts(camera_frame))
    front_end = ZernikeOpticalFrontEnd(spot, propagation, binning, nothing)
    wfs = ZernikeWFS{typeof(front_end),typeof(acquisition),typeof(selector)}(
        front_end, acquisition)
    initial_pupil = PupilFunction(tel)
    build_zernike_phasor!(zernike_propagation_workspace(wfs).phasor)
    build_zernike_phase_mask!(wfs, initial_pupil)
    return wfs
end

sensing_mode(::ZernikeWFS) = Diffractive()

@inline function _require_zernike_pupil_geometry(wfs::ZernikeWFS,
    pupil::PupilFunction)
    resolution = zernike_propagation_plan(wfs).pupil_resolution
    pupil.metadata.dimensions == (resolution, resolution) || throw(
        DimensionMismatchError(
            "ZernikeWFS PupilFunction dimensions do not match its prepared pupil grid"))
    return nothing
end

function build_zernike_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    return build_zernike_phasor!(execution_style(phasor), phasor)
end

function build_zernike_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for j in 1:n, i in 1:n
        phasor[i, j] = cis(scale * (i + j - 2))
    end
    return phasor
end

function build_zernike_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, zernike_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function host_zernike_phase_mask(wfs::ZernikeWFS, pupil::PupilFunction)
    phase_spot = wfs.front_end.phase_spot
    workspace = zernike_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(workspace.phase_mask, 1)
    T = eltype(zernike_acquisition_products(wfs).frame)
    host = Matrix{Complex{T}}(undef, pad, pad)
    center = T(pad) / 2
    radius = phase_spot.radius_lambda_over_d * (T(pad) / T(n))
    phase = T(pi) * phase_spot.phase_shift_pi
    shifted = cis(phase)
    @inbounds for j in 1:pad, i in 1:pad
        x = T(i) - center - T(0.5)
        y = T(j) - center - T(0.5)
        host[i, j] = hypot(x, y) <= radius ? shifted : one(shifted)
    end
    return host
end

function build_zernike_phase_mask!(wfs::ZernikeWFS, pupil::PupilFunction)
    phase_mask = zernike_propagation_workspace(wfs).phase_mask
    copyto!(phase_mask, host_zernike_phase_mask(wfs, pupil))
    return phase_mask
end

function sample_zernike_frame!(out::AbstractMatrix{T}, nominal::AbstractMatrix{T}, wfs::ZernikeWFS,
    input::AbstractMatrix{T}, pupil::PupilFunction) where {T<:AbstractFloat}
    plan = zernike_propagation_plan(wfs)
    sub = div(_pupil_resolution(pupil), plan.pupil_samples)
    bin2d!(nominal, input, sub)
    binning = zernike_acquisition_plan(wfs).binning
    if binning == 1
        copyto!(out, nominal)
    else
        bin2d!(out, nominal, binning)
    end
    return out
end

function zernike_pupil_intensity!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS")
    _require_zernike_pupil_geometry(wfs, pupil)
    propagation = zernike_propagation_workspace(wfs)
    T = eltype(zernike_acquisition_products(wfs).frame)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(
        photon_irradiance(src) * (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2
    ))
    amplitude = pupil.amplitude
    fill!(propagation.field, zero(eltype(propagation.field)))
    @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * amplitude *
        cispi(opd_to_cycles * pupil.opd)
    copyto!(propagation.focal_field, propagation.field)
    @. propagation.focal_field *= propagation.phasor
    execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
    @. propagation.focal_field *= propagation.phase_mask
    copyto!(propagation.pupil_field, propagation.focal_field)
    execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    @views @. propagation.pupil_intensity =
        abs2(propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    return propagation.pupil_intensity
end

@inline supports_prepared_runtime(::ZernikeWFS, src::AbstractSource) =
    is_leaf_source(src)
@inline supports_detector_output(::ZernikeWFS, ::AbstractDetector) = true

@inline function prepare_runtime_wfs!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS runtime preparation")
    _require_zernike_pupil_geometry(wfs, pupil)
    return wfs
end

include("zernike/stages.jl")

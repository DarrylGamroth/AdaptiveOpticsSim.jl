module ProperHILCoronagraphCommon

using AdaptiveOpticsSim
using Proper
using Random

export CoronagraphPayload,
    ProperHILCoronagraphContext,
    build_proper_hil_context,
    resolve_hil_backend,
    ao_step!,
    science_step!,
    hil_step!

struct CoronagraphPayload{P,A,T}
    opd_m::P
    pupil_amplitude::A
    diameter_m::T
    focal_length_m::T
    lyot_stop_norm::T
end

mutable struct ProperHILCoronagraphContext{
    TEL,S,A,R,P,O1,O2,W,PM,CP,SP,C1,C2,RNG,T,
}
    telescope::TEL
    source::S
    atmosphere::A
    atmosphere_renderer::R
    wfs_pupil::P
    tiptilt::O1
    dm::O2
    wfs::W
    science_model::PM
    payload::CP
    science_pupil::SP
    tiptilt_command::C1
    dm_command::C2
    rng::RNG
    step_index::Int
    wavelength_um::T
end

function resolve_hil_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return AdaptiveOpticsSim.CPUBackend()
    elseif lowered == "cuda"
        Base.find_package("CUDA") === nothing &&
            error("CUDA.jl is not available. Install it in the active environment before running the CUDA Proper HIL example.")
        @eval using CUDA
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        return AdaptiveOpticsSim.CUDABackend()
    elseif lowered == "amdgpu"
        Base.find_package("AMDGPU") === nothing &&
            error("AMDGPU.jl is not available. Install it in the active environment before running the AMDGPU Proper HIL example.")
        @eval using AMDGPU
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        return AdaptiveOpticsSim.AMDGPUBackend()
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function hil_coronagraph_prescription(λm, n; payload::CoronagraphPayload)
    wf = prop_begin(payload.diameter_m, λm, n; beam_diam_fraction=1.0)
    prop_circular_aperture(wf, payload.diameter_m / 2)
    prop_define_entrance(wf)
    prop_multiply(wf, payload.pupil_amplitude)
    prop_add_phase(wf, payload.opd_m)

    prop_lens(wf, payload.focal_length_m, "science arm")
    prop_propagate(wf, payload.focal_length_m, "focal plane mask")
    prop_8th_order_mask(wf, 4.0; circular=true)

    prop_propagate(wf, payload.focal_length_m, "lyot relay")
    prop_lens(wf, payload.focal_length_m, "lyot relay")
    prop_propagate(wf, 2payload.focal_length_m, "lyot stop")
    prop_circular_aperture(wf, payload.lyot_stop_norm; NORM=true)

    prop_propagate(wf, payload.focal_length_m, "science camera")
    prop_lens(wf, payload.focal_length_m, "science camera")
    prop_propagate(wf, payload.focal_length_m, "final focus")
    return prop_end(wf)
end

function _build_wfs(tel::Telescope; T::Type{<:AbstractFloat}, backend::AbstractArrayBackend)
    return ShackHartmannWFS(tel; n_lenslets=16, mode=Diffractive(), threshold=T(0), T=T, backend=backend)
end

function build_proper_hil_context(;
    backend::AbstractArrayBackend=AdaptiveOpticsSim.CPUBackend(),
    T::Type{<:AbstractFloat}=Float32,
    resolution::Int=128,
    diameter::Real=8.0,
    science_wavelength_um::Real=1.65,
    focal_length_m::Real=80.0,
    lyot_stop_norm::Real=0.9,
    rng_seed::Integer=1,
)
    selector = backend
    tel = Telescope(
        resolution=resolution,
        diameter=diameter,
        central_obstruction=T(0.1),
        T=T,
        backend=selector,
    )
    src = Source(band=:H, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0))
    tiptilt = TipTiltMirror(tel; scale=T(0.05), T=T, backend=selector, label=:tiptilt)
    dm = DeformableMirror(tel; n_act=16, influence_width=T(0.3), T=T, backend=selector)
    wfs = _build_wfs(tel; T=T, backend=selector)
    wfs_pupil = PupilFunction(tel)
    prepare_runtime_wfs!(wfs, wfs_pupil, src)
    atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, src)
    science_pupil = PupilFunction(tel)
    proper_ctx = Proper.RunContext(typeof(opd_map(science_pupil)))
    science_model = prepare_model(
        :proper_hil_coronagraph,
        hil_coronagraph_prescription,
        science_wavelength_um,
        resolution;
        context=proper_ctx,
        precision=T,
        pool_size=1,
    )

    payload = CoronagraphPayload(
        opd_map(science_pupil),
        pupil_amplitude(science_pupil),
        T(tel.params.diameter),
        T(focal_length_m),
        T(lyot_stop_norm),
    )

    tiptilt_command = similar(tiptilt.state.coefs)
    dm_command = similar(dm.state.coefs)
    fill!(tiptilt_command, zero(T))
    fill!(dm_command, zero(T))

    return ProperHILCoronagraphContext(
        tel,
        src,
        atm,
        atmosphere_renderer,
        wfs_pupil,
        tiptilt,
        dm,
        wfs,
        science_model,
        payload,
        science_pupil,
        tiptilt_command,
        dm_command,
        runtime_rng(rng_seed),
        0,
        T(science_wavelength_um),
    )
end

function _stage_command!(ctx::ProperHILCoronagraphContext)
    ctx.step_index += 1
    T = eltype(ctx.dm_command)
    tstep = T(ctx.step_index)
    ctx.tiptilt_command .= (
        T(5e-3) * sin(T(0.1) * tstep),
        -T(5e-3) * cos(T(0.1) * tstep),
    )
    fill!(ctx.dm_command, T(5e-9) * tstep)
    set_command!(ctx.tiptilt, ctx.tiptilt_command)
    set_command!(ctx.dm, ctx.dm_command)
    return (tiptilt=ctx.tiptilt_command, dm=ctx.dm_command)
end

function ao_step!(ctx::ProperHILCoronagraphContext)
    _stage_command!(ctx)
    epoch = advance_by!(ctx.atmosphere, 1e-3; rng=ctx.rng)
    render_atmosphere!(
        ctx.wfs_pupil,
        ctx.atmosphere_renderer,
        ctx.atmosphere,
        epoch,
    )
    for optic in (ctx.tiptilt, ctx.dm)
        update_surface!(optic)
        apply_surface!(ctx.wfs_pupil, optic, DMAdditive())
    end
    measure!(ctx.wfs, ctx.wfs_pupil, ctx.source)
    return slopes(ctx.wfs)
end

function science_step!(ctx::ProperHILCoronagraphContext)
    copyto!(ctx.science_pupil.amplitude, ctx.wfs_pupil.amplitude)
    copyto!(ctx.science_pupil.opd, ctx.wfs_pupil.opd)
    return prop_run(ctx.science_model; payload=ctx.payload)
end

function hil_step!(ctx::ProperHILCoronagraphContext)
    ao_step!(ctx)
    return science_step!(ctx)
end

end

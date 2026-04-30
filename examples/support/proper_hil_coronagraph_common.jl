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

struct CoronagraphPayload{P,M,T}
    phase_map_m::P
    pupil_mask::M
    diameter_m::T
    focal_length_m::T
    lyot_stop_norm::T
end

mutable struct ProperHILCoronagraphContext{SC,SM,PM,P,CB,R1,R2,T}
    scenario::SC
    sim::SM
    science_model::PM
    payload::P
    command_buffer::CB
    tiptilt_range::R1
    dm_range::R2
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

function hil_coronagraph_prescription(λm, n, payload::CoronagraphPayload)
    wf = prop_begin(payload.diameter_m, λm, n; beam_diam_fraction=1.0)
    prop_circular_aperture(wf, payload.diameter_m / 2)
    prop_define_entrance(wf)
    prop_multiply(wf, payload.pupil_mask)
    prop_add_phase(wf, payload.phase_map_m)

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

function _segment_range(layout::AdaptiveOpticsSim.RuntimeCommandLayout, label::Symbol)
    for seg in command_segments(layout)
        seg.label == label && return command_segment_range(seg)
    end
    throw(InvalidConfiguration("missing command segment '$label' in runtime command layout"))
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
        sampling_time=T(1e-3),
        central_obstruction=T(0.1),
        T=T,
        backend=selector,
    )
    src = Source(band=:H, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0))
    tiptilt = TipTiltMirror(tel; scale=T(0.05), T=T, backend=selector, label=:tiptilt)
    dm = DeformableMirror(tel; n_act=16, influence_width=T(0.3), T=T, backend=selector)
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = _build_wfs(tel; T=T, backend=selector)
    sim = AOSimulation(tel, src, atm, optic, wfs)

    branch = ControlLoopBranch(:main, sim, NullReconstructor(); rng=runtime_rng(rng_seed))
    cfg = SingleControlLoopConfig(
        name=:proper_hil,
        branch_label=:main,
        outputs=RuntimeOutputRequirements(slopes=true),
    )
    scenario = build_control_loop_scenario(cfg, branch)
    prepare!(scenario)

    proper_ctx = Proper.RunContext(typeof(sim.tel.state.opd))
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
        sim.tel.state.opd,
        sim.tel.state.pupil,
        T(sim.tel.params.diameter),
        T(focal_length_m),
        T(lyot_stop_norm),
    )

    command_buffer = similar(command(scenario))
    copyto!(command_buffer, command(scenario))
    layout = AdaptiveOpticsSim.command_layout(scenario)
    tiptilt_range = _segment_range(layout, :tiptilt)
    dm_range = _segment_range(layout, :dm)

    return ProperHILCoronagraphContext(
        scenario,
        sim,
        science_model,
        payload,
        command_buffer,
        tiptilt_range,
        dm_range,
        0,
        T(science_wavelength_um),
    )
end

function _stage_command!(ctx::ProperHILCoronagraphContext)
    ctx.step_index += 1
    T = eltype(ctx.command_buffer)
    phase = T(ctx.step_index)
    @views begin
        tiptilt = ctx.command_buffer[ctx.tiptilt_range]
        dm = ctx.command_buffer[ctx.dm_range]
        tiptilt[1] = T(5e-3) * sin(T(0.1) * phase)
        tiptilt[2] = -T(5e-3) * cos(T(0.1) * phase)
        fill!(dm, T(5e-9) * phase)
    end
    set_command!(ctx.scenario, ctx.command_buffer)
    return ctx.command_buffer
end

function ao_step!(ctx::ProperHILCoronagraphContext)
    _stage_command!(ctx)
    sense!(ctx.scenario)
    return slopes(ctx.scenario)
end

function science_step!(ctx::ProperHILCoronagraphContext)
    return prop_run(ctx.science_model; PASSVALUE=ctx.payload)
end

function hil_step!(ctx::ProperHILCoronagraphContext)
    ao_step!(ctx)
    return science_step!(ctx)
end

end

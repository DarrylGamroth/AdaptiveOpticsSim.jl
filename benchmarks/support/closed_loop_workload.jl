mutable struct ClosedLoopWorkload{
    S,A,R,P,O,W,C,V,RNG,T,
}
    source::S
    atmosphere::A
    atmosphere_renderer::R
    pupil::P
    optic::O
    wfs::W
    reconstructor::C
    command::V
    rng::RNG
    sample_period::T
end

function prepare_closed_loop_workload(;
    resolution::Int=16,
    n_lenslets::Int=4,
    n_act::Int=4,
    sample_period::Real=1e-3,
    T::Type{<:AbstractFloat}=Float32,
    backend::AbstractArrayBackend=CPUBackend(),
    seed::Integer=1,
)
    telescope = Telescope(
        resolution=resolution,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    source = Source(band=:I, magnitude=0, T=T)
    atmosphere = KolmogorovAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        T=T,
        backend=backend,
    )
    optic = DeformableMirror(
        telescope;
        n_act=n_act,
        influence_width=T(0.3),
        T=T,
        backend=backend,
    )
    wfs = ShackHartmannWFS(
        telescope;
        n_lenslets=n_lenslets,
        mode=Diffractive(),
        T=T,
        backend=backend,
    )
    source_pupil = PupilFunction(telescope; T=T, backend=backend)
    interaction = interaction_matrix(
        optic,
        wfs,
        source_pupil,
        source;
        amplitude=T(0.05),
    )
    reconstructor = ModalReconstructor(interaction; gain=T(0.5))
    command = similar(optic.state.coefs)
    fill!(command, zero(T))
    atmosphere_renderer = prepare_atmosphere_renderer(
        atmosphere,
        telescope,
        source,
    )
    workload = ClosedLoopWorkload(
        source,
        atmosphere,
        atmosphere_renderer,
        source_pupil,
        optic,
        wfs,
        reconstructor,
        command,
        runtime_rng(seed),
        T(sample_period),
    )
    step_closed_loop_workload!(workload)
    return workload
end

function step_closed_loop_workload!(workload::ClosedLoopWorkload)
    epoch = advance_by!(
        workload.atmosphere,
        workload.sample_period;
        rng=workload.rng,
    )
    render_atmosphere!(
        workload.pupil,
        workload.atmosphere_renderer,
        workload.atmosphere,
        epoch,
    )
    update_surface!(workload.optic)
    apply_surface!(workload.pupil, workload.optic, DMAdditive())
    measure!(workload.wfs, workload.pupil, workload.source)
    reconstruct!(
        workload.command,
        workload.reconstructor,
        slopes(workload.wfs),
    )
    @. workload.command = -workload.command
    copyto!(workload.optic.state.coefs, workload.command)
    return workload
end

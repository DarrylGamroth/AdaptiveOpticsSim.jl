function bind_camera_opd(pupil::O.PupilFunction{M,S,A,I,Backend}, opd) where {M,S,A,I,Backend}
    Base.require_one_based_indexing(opd)
    O.validate_plane_storage(pupil.metadata, opd; label="camera OPD")
    return O.PupilFunction{M,S,A,typeof(opd),Backend}(
        pupil.metadata, pupil.support, pupil.amplitude, opd, pupil.aperture_revision,
    )
end

function FG.prepare_cpu_loan_node(prepared::FGA.PreparedAlgorithm{PyramidRateF32},
        outputs::Tuple, inputs::Tuple, ::Val)
    opd, photon_rate = only(inputs), only(outputs)
    owner = prepared.workspace.owner
    optics = owner.prepared
    pupil = bind_camera_opd(optics.input, opd)
    output = O.IntensityMap(optics.output.metadata, photon_rate)
    # Only the pupil and photon-rate bindings differ across slots. Reuse the
    # same propagation plan, FFT plans, and workspace under serialized replay.
    bound = AdaptiveOpticsSim.WavefrontSensors.PreparedPyramidOptics(
        optics.plan, optics.front_end, optics.workspace, pupil, output,
        optics.workspace_binding, optics.backend, optics.device,
    )
    bound_owner = AOG._PyramidRateOwner(bound, owner.field_input, owner.amplitude_scale)
    workspace = PyramidRateWorkspace(bound_owner, opd, photon_rate)
    return FGA.PreparedAlgorithm(PyramidRateF32, prepared.plan, workspace,
        FGA.algorithm_formats(prepared))
end

function FG.prepare_cpu_loan_node(prepared::FGA.PreparedAlgorithm{EMCCDF32},
        outputs::Tuple, inputs::Tuple, ::Val{true})
    # Internal ports retain their addresses across every slot, so one detector
    # owns the seeded stream for the whole graph.
    bound = prepare_host_emccd(prepared.plan, only(inputs), only(outputs))
    return FGA.PreparedAlgorithm(EMCCDF32, bound.plan, bound.workspace,
        FGA.algorithm_formats(prepared))
end

function bind_shack_hartmann_loan(optics::AdaptiveOpticsSim.WavefrontSensors.PreparedShackHartmannOptics,
        opd, photon_rate)
    pupil = bind_camera_opd(optics.input, opd)
    output = O.IntensityMap(optics.output.metadata, photon_rate)
    return AdaptiveOpticsSim.WavefrontSensors.PreparedShackHartmannOptics(
        optics.plan, optics.optics, optics.workspace, pupil, output,
        optics.workspace_binding, optics.layout_binding, optics.backend, optics.device,
    )
end

# Spectral bundles retain the existing copy-based binding until their component
# inputs have been qualified for rotating external loans.
bind_shack_hartmann_loan(optics, opd, photon_rate) = nothing

function FG.prepare_cpu_loan_node(prepared::FGA.PreparedAlgorithm{ShackHartmannRateF32},
        outputs::Tuple, inputs::Tuple, ::Val)
    opd, photon_rate = only(inputs), only(outputs)
    optics = bind_shack_hartmann_loan(prepared.workspace.owner.prepared, opd, photon_rate)
    optics === nothing && return prepared
    workspace = ShackHartmannRateWorkspace(AOG._ShackHartmannRateOwner(optics), opd, photon_rate)
    return FGA.PreparedAlgorithm(ShackHartmannRateF32, prepared.plan, workspace,
        FGA.algorithm_formats(prepared))
end

function FG.prepare_cpu_loan_node(prepared::FGA.PreparedAlgorithm{CMOSF32},
        outputs::Tuple, inputs::Tuple, ::Val{true})
    bound = prepare_host_cmos(prepared.plan, only(inputs), only(outputs))
    return FGA.PreparedAlgorithm(CMOSF32, bound.plan, bound.workspace,
        FGA.algorithm_formats(prepared))
end

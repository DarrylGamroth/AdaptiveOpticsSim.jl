include(joinpath(@__DIR__, "common.jl"))

using AdaptiveOpticsCalibration

const PR = AdaptiveOpticsCalibration.PhaseRetrieval

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    zb = ZernikeBasis(tel, 6)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 1:4]
    coeffs_true = [25e-9, -10e-9, 5e-9, 0.0]
    pupil = PupilFunction(tel)
    apply_opd!(pupil, combine_modes(basis, coeffs_true))
    imaging = prepare_direct_imaging(pupil, src; zero_padding=zero_padding)
    image = copy(intensity_values(form_direct_image!(imaging)))
    diversity = similar(pupil.opd)
    fill!(diversity, zero(eltype(diversity)))
    forward = prepare_lift_forward_model(tel, src, basis, pupil.opd;
        diversity_opd=diversity, zero_padding=zero_padding)
    observation = LiFTObservation(forward, image)
    specification = PR.LiFTSpecification(forward, observation)
    plan = AdaptiveOpticsCalibration.prepare(PR.LiFT(iterations=3,
            jacobian_method=PR.LiFTAnalyticJacobian()), specification)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = PR.LiFTInputs(observation.values)
    AdaptiveOpticsCalibration.process!(result, workspace, plan, inputs)
    coeffs_fit = PR.lift_coefficients(result)

    @info "LiFT tutorial complete" n_modes=length(coeffs_true)
    return (
        coeffs_true=coeffs_true,
        coeffs_fit=coeffs_fit,
        image=image,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

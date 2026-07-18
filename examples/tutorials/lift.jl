include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    zb = ZernikeBasis(tel, 6)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 1:4]
    coeffs_true = [25e-9, -10e-9, 5e-9, 0.0]
    apply_opd!(tel, combine_modes(basis, coeffs_true))

    pupil = PupilFunction(tel)
    apply_opd!(pupil, opd_map(tel))
    imaging = prepare_direct_imaging(tel, pupil, src;
        zero_padding=zero_padding)
    image = copy(intensity_values(form_direct_image!(imaging)))
    diversity = zeros(eltype(tel.state.opd), size(tel.state.opd))
    forward = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, zero_padding=zero_padding)
    lift = LiFT(forward; iterations=3, numerical=false)
    observation = LiFTObservation(forward, image)
    coeffs_fit = AdaptiveOpticsSim.reconstruct(lift, observation;
        coeffs0=zeros(4))

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

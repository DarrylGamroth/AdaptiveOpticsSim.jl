include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    zb = ZernikeBasis(tel, 6)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 1:4]
    coeffs_true = [25e-9, -10e-9, 5e-9, 0.0]
    apply_opd!(tel, combine_modes(basis, coeffs_true))

    psf = copy(compute_psf!(tel, src; zero_padding=zero_padding))
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=zero_padding, binning=1)
    diversity = zeros(eltype(tel.state.opd), size(tel.state.opd))
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=3, numerical=false)
    coeffs_fit = AdaptiveOpticsSim.reconstruct(lift, psf, collect(1:length(coeffs_true)); coeffs0=zeros(4))

    @info "LiFT tutorial complete" n_modes=length(coeffs_true)
    return (
        coeffs_true=coeffs_true,
        coeffs_fit=coeffs_fit,
        psf=psf,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

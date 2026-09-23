include(joinpath(@__DIR__, "common.jl"))
import AdaptiveOpticsCalibration

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    coefficients = [0.0, 30e-9, -20e-9, 10e-9]
    zernike = ZernikeBasis(tel, length(coefficients))
    compute_zernike!(zernike, tel)
    modal = AdaptiveOpticsCalibration.ModalBases
    specification = modal.ModalOPDExpansionSpecification(
        resolution, resolution, length(coefficients), pupil_mask(tel), Float64)
    plan = AdaptiveOpticsCalibration.prepare(modal.ModalOPDExpansion(), specification)
    product = AdaptiveOpticsCalibration.process(
        plan, modal.ModalOPDExpansionInputs(zernike.modes, coefficients))
    ncpa = NCPA(product.opd)
    pupil = PupilFunction(tel)
    apply_surface!(pupil, ncpa, DMReplace())
    imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
    image = copy(intensity_values(form_direct_image!(imaging)))

    @info "NCPA tutorial complete" opd_rms=pupil_rms(pupil.opd, pupil_support(pupil))
    return (
        ncpa_opd=copy(pupil.opd),
        image=image,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

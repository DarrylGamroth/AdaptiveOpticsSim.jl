include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(magnitude=6.0)
    model = GaussianDiskSourceModel(sigma_arcsec=0.4, n_side=5)
    ext = with_extended_source(src, model)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @. tel.state.opd = 5e-8 * zb.modes[:, :, 5]

    sh_point = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive())
    sh_ext = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive())
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_point, tel, src)
    ext_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext, tel, ext)
    point_slopes = copy(measure!(sh_point, tel, src))
    ext_slopes = copy(measure!(sh_ext, tel, ext))

    pyr_point = PyramidWFS(tel; pupil_samples=6, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=6, mode=Diffractive(), modulation=1.0)
    pyr_point_slopes = copy(measure!(pyr_point, tel, src))
    pyr_ext_slopes = copy(measure!(pyr_ext, tel, ext))

    @info "Extended-source sensing tutorial complete" sh_delta_norm=LinearAlgebra.norm(sh_ext.state.spot_cube .- sh_point.state.spot_cube) pyramid_delta_norm=LinearAlgebra.norm(pyr_ext.state.intensity .- pyr_point.state.intensity)
    return (
        sh_point_peak=point_peak,
        sh_extended_peak=ext_peak,
        sh_point_slopes=point_slopes,
        sh_extended_slopes=ext_slopes,
        sh_spot_delta=copy(sh_ext.state.spot_cube .- sh_point.state.spot_cube),
        pyramid_point_slopes=pyr_point_slopes,
        pyramid_extended_slopes=pyr_ext_slopes,
        pyramid_intensity_delta=copy(pyr_ext.state.intensity .- pyr_point.state.intensity),
        n_samples=length(extended_source_asterism(ext)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

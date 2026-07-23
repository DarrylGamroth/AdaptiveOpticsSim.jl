using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

# This standalone audit intentionally exercises package-internal backend seams.
# Keep those names local to the script without expanding the public API.
for name in names(AdaptiveOpticsSim; all=true)
    text = String(name)
    if Base.isidentifier(text) && !startswith(text, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

function record_gpu_smoke!(f::Function, failures::Vector{String}, name::AbstractString)
    try
        result = f()
        println(name, " ok :: ", typeof(result))
    catch err
        push!(failures, string(name, " ERROR :: ", sprint(showerror, err)))
        println(name, " ERROR :: ", sprint(showerror, err))
    end
    return nothing
end

function prepared_gpu_field(tel::Telescope, src::AbstractSource;
    zero_padding::Int, T::Type{<:AbstractFloat})
    wavefront = PupilFunction(tel; T=T)
    field = ElectricField(wavefront, src; zero_padding=zero_padding, T=T)
    plan = prepare_pupil_field(wavefront, src, field)
    fill_electric_field!(field, wavefront, plan)
    return (; wavefront, field, plan)
end

function gpu_direct_image(tel::Telescope, src::AbstractSource;
    zero_padding::Int, T::Type{<:AbstractFloat})
    pupil = PupilFunction(tel; T=T)
    prepared = prepare_direct_imaging(pupil, src; zero_padding=zero_padding)
    return form_direct_image!(prepared)
end

function run_gpu_smoke_matrix(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    failures = String[]
    rng = MersenneTwister(1)
    T = Float32
    atmosphere_step = T(1e-3)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    backend = AdaptiveOpticsSim.array_backend_selector(BackendArray)

    tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    lgs = LGSSource(; magnitude=0.0, wavelength=589e-9, altitude=90_000.0,
        laser_coordinates=(0.0, 0.0), photon_irradiance=one(T), T=T)
    spider_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
    apply_spiders!(spider_tel; thickness=0.5, angles=[0.0, 90.0])
    pupil = PupilFunction(tel; T=T, backend=backend)
    spider_pupil = PupilFunction(spider_tel; T=T, backend=backend)

    record_gpu_smoke!(failures, "direct_image_source") do
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        @assert rate_map.values isa BackendArray
        return rate_map.values
    end

    record_gpu_smoke!(failures, "direct_image_off_axis_orientation") do
        image_pupil = PupilFunction(tel; T=T)
        on_axis = prepare_direct_imaging(image_pupil, src; zero_padding=2)
        sample_arcsec = focal_plane_pixel_scale_arcsec(
            direct_imaging_output(on_axis))
        positive_x = Source(band=:I, magnitude=zero(T),
            coordinates=(T(sample_arcsec), zero(T)), T=T)
        off_axis = prepare_direct_imaging(image_pupil, positive_x;
            zero_padding=2)
        @assert off_axis.plan.shift_samples == (1, 0)
        reference = Array(intensity_values(form_direct_image!(on_axis)))
        shifted = Array(intensity_values(form_direct_image!(off_axis)))
        relative_error = maximum(abs.(shifted .-
            circshift(reference, (1, 0)))) /
            max(maximum(abs, reference), eps(Float64))
        @assert relative_error < 2f-5
        return off_axis.output.values
    end

    record_gpu_smoke!(failures, "electric_field_core") do
        prepared = prepared_gpu_field(tel, src; zero_padding=2, T=T)
        field = prepared.field
        @assert field.values isa BackendArray
        intensity = similar(field.values, T)
        intensity!(intensity, field)
        @assert intensity isa BackendArray

        amplitude = similar(prepared.wavefront.opd, T,
            tel.params.resolution, tel.params.resolution)
        fill!(amplitude, T(0.5))
        apply_amplitude!(field, amplitude, prepared.plan)
        intensity!(intensity, field)
        @assert intensity isa BackendArray

        prepared2 = prepared_gpu_field(tel, src; zero_padding=2, T=T)
        field2 = prepared2.field
        direct_from_field = prepare_direct_imaging(src, field2)
        rate_from_field = form_direct_image!(direct_from_field)
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        rel = maximum(abs.(Array(rate_from_field.values) .-
            Array(rate_map.values))) /
            max(maximum(abs.(Array(rate_map.values))), eps(Float64))
        @assert rel < 1f-5
        return intensity
    end

    record_gpu_smoke!(failures, "field_propagation") do
        prepared = prepared_gpu_field(tel, src; zero_padding=2, T=T)
        field = prepared.field
        fraunhofer = FraunhoferPropagation(field)
        propagated = propagation_output(field, fraunhofer)
        propagate_field!(propagated, field, fraunhofer)
        rate_from_propagation = similar(field.values, T)
        @. rate_from_propagation = abs2(propagated.values)
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        rel = maximum(abs.(Array(rate_from_propagation) .-
            Array(rate_map.values))) /
            max(maximum(abs.(Array(rate_map.values))), eps(Float64))
        @assert rel < 1f-5

        fresnel = FresnelPropagation(field; distance_m=T(10))
        fresnel_out = propagation_output(field, fresnel)
        propagate_field!(fresnel_out, field, fresnel)
        @assert fresnel_out.values isa BackendArray
        reverse = FresnelPropagation(fresnel_out; distance_m=T(-10),
            output_kind=PupilPlane())
        field_reverse = propagation_output(fresnel_out, reverse)
        propagate_field!(field_reverse, fresnel_out, reverse)
        roundtrip_rel = maximum(abs.(Array(field_reverse.values) .-
            Array(field.values))) /
            max(maximum(abs.(Array(field.values))), eps(Float64))
        @assert roundtrip_rel < 1f-4
        return fresnel_out.values
    end

    record_gpu_smoke!(failures, "direct_image_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0)),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0)),
        ])
        rate_map = gpu_direct_image(tel, ast; zero_padding=2, T=T)
        @assert rate_map.values isa BackendArray
        return rate_map.values
    end

    record_gpu_smoke!(failures, "direct_image_source_spiders") do
        rate_map = gpu_direct_image(spider_tel, src; zero_padding=2, T=T)
        @assert rate_map.values isa BackendArray
        return rate_map.values
    end

    record_gpu_smoke!(failures, "aperture_masks") do
        support = BackendArray{Bool}(undef, 16, 16)
        weighted = BackendArray{Complex{T}}(undef, 16, 16)
        valid = BackendArray{Bool}(undef, 4, 4)
        build_mask!(support, AnnularAperture(inner_radius=T(0.2), outer_radius=T(1), T=T); grid=default_mask_grid(support; T=T))
        apply_mask!(support, SpiderMask(thickness=T(0.1), angle_rad=T(pi / 4), T=T); grid=default_mask_grid(support; T=T))
        build_mask!(weighted, CircularAperture(radius=T(4), T=T); grid=pixel_mask_grid(weighted; T=T), inside=complex(T(inv(sqrt(2))), T(inv(sqrt(2)))))
        build_mask!(valid, SubapertureGridMask(threshold=T(0.1), T=T), support)
        @assert support isa BackendArray
        @assert weighted isa BackendArray
        @assert valid isa BackendArray
        @assert any(Array(valid))
        return valid
    end

    record_gpu_smoke!(failures, "detector_capture_none") do
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=2, T=T, backend=backend)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(det, rate_map, acquisition; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "detector_capture_noise") do
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        det = Detector(noise=(NoisePhoton(), NoiseReadout(T(1e-3))), integration_time=1.0, qe=1.0,
            binning=2, background_flux=T(0.5), dark_current=T(0.1), T=T, backend=backend)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(det, rate_map, acquisition; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "atmosphere_step") do
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=backend)
        output = PupilFunction(tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, tel)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_multilayer_step") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        output = PupilFunction(tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, tel)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert atm.layers[1].generator.state.opd isa BackendArray
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_multilayer_step_spiders") do
        atm = MultiLayerAtmosphere(spider_tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        output = PupilFunction(spider_tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, spider_tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert atm.layers[1].generator.state.opd isa BackendArray
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_multilayer_step") do
        atm = InfiniteMultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
            T=T,
            backend=backend,
        )
        output = PupilFunction(tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert output.opd isa BackendArray
        @assert atm.layers[1].screen.state.screen isa BackendArray
        @assert atm.layers[1].screen.state.screen_scratch isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_multilayer_step_spiders") do
        atm = InfiniteMultiLayerAtmosphere(spider_tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
            T=T,
            backend=backend,
        )
        output = PupilFunction(spider_tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, spider_tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmospheric_field_geometric") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        field_pupil = PupilFunction(tel; T=T, backend=backend)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        prop = AtmosphericFieldPropagation(atm, field_pupil, src;
            model=GeometricAtmosphericPropagation(T=T),
            zero_padding=2,
            T=T)
        field = propagate_atmosphere_field!(prop, atm, epoch)
        @assert field.values isa BackendArray
        intensity = atmospheric_intensity!(prop, atm, epoch)
        @assert intensity isa BackendArray
        return intensity
    end

    record_gpu_smoke!(failures, "atmospheric_field_fresnel") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        field_pupil = PupilFunction(tel; T=T, backend=backend)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        prop = AtmosphericFieldPropagation(atm, field_pupil, src;
            model=LayeredFresnelAtmosphericPropagation(T=T),
            zero_padding=2,
            T=T)
        field = propagate_atmosphere_field!(prop, atm, epoch)
        @assert field.values isa BackendArray
        return field.values
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_statistical_agreement") do
        function trajectory_stats(backend; steps::Int=10)
            local_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
            local_src = Source(band=:I, magnitude=0.0, coordinates=(30.0, 20.0), T=T)
            local_atm = InfiniteMultiLayerAtmosphere(local_tel;
                r0=T(0.2),
                L0=T(25.0),
                fractional_cn2=T[0.7, 0.3],
                wind_speed=T[8.0, 4.0],
                wind_direction=T[0.0, 90.0],
                altitude=T[0.0, 5000.0],
                screen_resolution=33,
                stencil_size=35,
                T=T,
                backend=backend,
            )
            local_pupil = PupilFunction(local_tel; T=T, backend=backend)
            renderer = prepare_atmosphere_renderer(local_atm, local_tel,
                local_src)
            local_rng = MersenneTwister(21)
            stds = Float64[]
            corrs = Float64[]
            previous = nothing
            for _ in 1:steps
                epoch = advance_by!(local_atm, atmosphere_step; rng=local_rng)
                render_atmosphere!(local_pupil, renderer, local_atm, epoch)
                opd = Array(local_pupil.opd)
                push!(stds, std(vec(opd)))
                if previous !== nothing
                    push!(corrs, dot(vec(previous), vec(opd)) / sqrt(dot(vec(previous), vec(previous)) * dot(vec(opd), vec(opd))))
                end
                previous = opd
            end
            return (; std_mean=mean(stds), corr_mean=mean(corrs))
        end

        cpu_stats = trajectory_stats(CPUBackend())
        gpu_stats = trajectory_stats(backend)
        std_rel = abs(gpu_stats.std_mean - cpu_stats.std_mean) / max(cpu_stats.std_mean, eps(Float64))
        corr_abs = abs(gpu_stats.corr_mean - cpu_stats.corr_mean)
        @assert std_rel < 0.5
        @assert corr_abs < 0.05
        return gpu_stats
    end

    record_gpu_smoke!(failures, "atmosphere_phase_helpers") do
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=backend)
        ws = PhaseStatsWorkspace(tel.params.resolution; T=T, backend=backend)
        screen, psd = ft_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, return_psd=true)
        sh_screen = ft_sh_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, subharmonics=true, n_levels=2, subharmonic_radius=1)
        rho = similar(atm.state.opd, T, 4, 4)
        copyto!(rho, reshape(T[0.0, 0.02, 0.05, 0.1, 0.15, 0.25, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0], 4, 4))
        cov = phase_covariance(rho, atm)
        freqs = similar(atm.state.freqs, T, 4)
        copyto!(freqs, T[0.1, 0.2, 0.3, 0.4])
        covmat = covariance_matrix(freqs, freqs, atm)
        spectrum = phase_spectrum(freqs, atm)
        freq_grid = similar(atm.state.opd, T, 2, 2)
        copyto!(freq_grid, reshape(T[0.1, 0.2, 0.3, 0.4], 2, 2))
        spectrum_grid = phase_spectrum(freq_grid, atm)
        @assert screen isa BackendArray
        @assert psd isa BackendArray
        @assert sh_screen isa BackendArray
        @assert cov isa BackendArray
        @assert covmat isa BackendArray
        @assert spectrum isa BackendArray
        @assert spectrum_grid isa BackendArray
        cov_ref = phase_covariance(Array(rho), atm)
        covmat_ref = covariance_matrix(Array(freqs), Array(freqs), atm)
        cov_rel = maximum(abs.(Array(cov) .- cov_ref)) / maximum(abs.(cov_ref))
        covmat_rel = maximum(abs.(Array(covmat) .- covmat_ref)) / maximum(abs.(covmat_ref))
        @assert cov_rel < 5e-4
        @assert covmat_rel < 5e-4
        return screen
    end

    record_gpu_smoke!(failures, "dm_apply") do
        dm = DeformableMirror(tel; n_act=4, influence_width=0.3, T=T, backend=backend)
        dm_pupil = PupilFunction(tel; T=T, backend=backend)
        fill!(dm.state.coefs, T(0.05))
        update_surface!(dm)
        apply_surface!(dm_pupil, dm, DMReplace())
        @assert dm_pupil.opd isa BackendArray
        return dm_pupil.opd
    end

    record_gpu_smoke!(failures, "measure_shack_geometric") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        slopes = measure!(wfs, pupil)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_spectral_common_grid") do
        bundle = SpectralBundle(fill(wavelength(src), 2), T[0.4, 0.6]; T=T)
        spectral = with_spectrum(src, bundle)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, spectral)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "reject_shack_diffractive_distinct_wavelength_grids") do
        bundle = SpectralBundle(
            T[0.9 * wavelength(src), 1.1 * wavelength(src)],
            T[0.4, 0.6]; T=T)
        spectral = with_spectrum(src, bundle)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(),
            T=T, backend=backend)
        rejected = false
        try
            measure!(wfs, pupil, spectral)
        catch err
            err isa InvalidConfiguration || rethrow()
            rejected = true
        end
        @assert rejected
        return nothing
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_extended") do
        model = GaussianDiskSourceModel(sigma_arcsec=T(0.35), n_side=5, T=T)
        ext = with_extended_source(src, model)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, ext)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_spiders") do
        wfs = ShackHartmannWFS(spider_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, spider_pupil, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_lgs") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, lgs)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, ast)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_asterism_detector") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=backend)
        slopes = measure!(wfs, pupil, ast, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_detector_equivalence") do
        cpu_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=CPUBackend())
        gpu_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        cpu_src = Source(band=:I, magnitude=0.0, T=T)
        gpu_src = Source(band=:I, magnitude=0.0, T=T)
        cpu_wfs = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
        gpu_wfs = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        cpu_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
            sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
        gpu_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
            sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)

        cpu_pupil = PupilFunction(cpu_tel; T=T, backend=CPUBackend())
        gpu_pupil = PupilFunction(gpu_tel; T=T, backend=backend)
        measure!(cpu_wfs, cpu_pupil, cpu_src, cpu_det; rng=rng)
        measure!(gpu_wfs, gpu_pupil, gpu_src, gpu_det; rng=rng)

        cpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(cpu_wfs))
        gpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(gpu_wfs))
        cpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(cpu_wfs, cpu_det))
        gpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(gpu_wfs, gpu_det))

        @assert size(gpu_export) == size(cpu_export)
        @assert size(gpu_frame) == size(cpu_frame)
        @assert isapprox(gpu_export, cpu_export; rtol=1f-5, atol=1f-4)
        @assert isapprox(gpu_frame, cpu_frame; rtol=1f-5, atol=1f-4)
        return gpu_frame
    end

    record_gpu_smoke!(failures, "measure_pyramid_geometric") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T, backend=backend)
        slopes = measure!(wfs, pupil)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_diffractive") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_diffractive_polychromatic") do
        bundle = SpectralBundle(T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6]; T=T)
        poly = with_spectrum(src, bundle)
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, poly)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_diffractive_extended") do
        model = GaussianDiskSourceModel(sigma_arcsec=T(0.35), n_side=5, T=T)
        ext = with_extended_source(src, model)
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, ext)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_detector") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, mode=Diffractive(), T=T, backend=backend)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=backend)
        slopes = measure!(wfs, pupil, src, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_bioedge_geometric") do
        wfs = BioEdgeWFS(tel; pupil_samples=4, modulation=0.0, T=T, backend=backend)
        slopes = measure!(wfs, pupil)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_bioedge_diffractive") do
        wfs = BioEdgeWFS(tel; pupil_samples=4, modulation=2.0, mode=Diffractive(), T=T, backend=backend)
        slopes = measure!(wfs, pupil, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_zernike_diffractive") do
        wfs = ZernikeWFS(tel; pupil_samples=4, T=T, backend=backend)
        slopes = measure!(wfs, pupil, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_zernike_detector") do
        wfs = ZernikeWFS(tel; pupil_samples=4, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=backend)
        slopes = measure!(wfs, pupil, src, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_curvature_atmosphere") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        advance_by!(atm, atmosphere_step; rng=rng)
        wfs = CurvatureWFS(tel; pupil_samples=4, T=T, backend=backend)
        slopes = measure!(wfs, pupil, src, atm)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "closed_loop_step") do
        step_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        atm = KolmogorovAtmosphere(step_tel; r0=0.2, L0=25.0, T=T, backend=backend)
        dm = DeformableMirror(step_tel; n_act=4, influence_width=0.3, T=T, backend=backend)
        wfs = ShackHartmannWFS(step_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=backend)
        step_pupil = PupilFunction(step_tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, step_tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(step_pupil, renderer, atm, epoch)
        update_surface!(dm)
        apply_surface!(step_pupil, dm, DMAdditive())
        slopes = measure!(wfs, step_pupil, src, det; rng=rng)
        imaging = prepare_direct_imaging(step_pupil, src; zero_padding=2)
        rate_map = form_direct_image!(imaging)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(det, rate_map, acquisition; rng=rng)
        @assert slopes isa BackendArray
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "interaction_matrix_reconstructor") do
        cal_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        dm = DeformableMirror(cal_tel; n_act=4, influence_width=0.3, T=T, backend=backend)
        wfs = ShackHartmannWFS(cal_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
        calibration_pupil = PupilFunction(cal_tel; T=T, backend=backend)
        imat = interaction_matrix(dm, wfs, calibration_pupil, src;
            amplitude=T(0.05))
        recon = ModalReconstructor(imat; gain=one(T))
        measure!(wfs, calibration_pupil, src)
        cmds = reconstruct(recon, slopes(wfs))
        @assert imat.matrix isa BackendArray
        @assert cmds isa BackendArray
        return cmds
    end

    record_gpu_smoke!(failures, "gain_sensing_camera") do
        mask = backend_fill(B, one(T), 8, 8)
        basis = backend_rand(B, T, 8, 8, 3)
        frame = abs.(backend_randn(B, T, 8, 8))
        gsc = GainSensingCamera(mask, basis; T=T)
        calibrate!(gsc, frame)
        og = compute_optical_gains!(gsc, frame)
        @assert og isa BackendArray
        return og
    end

    record_gpu_smoke!(failures, "lift") do
        lift_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        lift_src = Source(band=:I, magnitude=8.0, T=T)
        basis = backend_rand(B, T, 16, 16, 3)
        diversity = backend_zeros(B, T, 16, 16)
        forward = prepare_lift_forward_model(lift_tel, lift_src, basis;
            diversity_opd=diversity, focal_resolution=32,
            zero_padding=2)
        lift = LiFT(forward; iterations=2, mode_ids=(1, 2),
            solve_mode=LiFTSolveAuto())
        rate_map = gpu_direct_image(lift_tel, lift_src;
            zero_padding=2, T=T)
        observation = LiFTObservation(forward, rate_map.values)
        coeffs = reconstruct(lift, observation)
        @assert coeffs isa BackendArray
        return coeffs
    end

    if !isempty(failures)
        error("GPU smoke matrix failed:\n" * join(failures, "\n"))
    end

    println("gpu_smoke_matrix complete")
    return nothing
end

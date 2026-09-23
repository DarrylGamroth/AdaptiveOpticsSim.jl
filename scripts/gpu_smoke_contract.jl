using AdaptiveOpticsSim
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsCalibration
using KernelAbstractions
using LinearAlgebra
using Random
using Statistics

const PR = AdaptiveOpticsCalibration.PhaseRetrieval

# This standalone audit intentionally exercises package-internal backend seams.
# Keep those names local to the script without expanding the public API.
for name in names(AdaptiveOpticsSim; all=true)
    text = String(name)
    if Base.isidentifier(text) && !startswith(text, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

for name in names(Backends; all=true)
    text = String(name)
    if Base.isidentifier(text) && !startswith(text, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Backends, $(QuoteNode(name)))
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

function prepared_gpu_field(tel::Telescope, src;
    zero_padding::Int, T::Type{<:AbstractFloat})
    wavefront = PupilFunction(tel; T=T)
    field = ElectricField(wavefront, src; zero_padding=zero_padding, T=T)
    plan = prepare_pupil_field(wavefront, src, field)
    fill_electric_field!(field, wavefront, plan)
    return (; wavefront, field, plan)
end

function gpu_direct_image(tel::Telescope, src;
    zero_padding::Int, T::Type{<:AbstractFloat})
    pupil = PupilFunction(tel; T=T)
    prepared = prepare_direct_imaging(pupil, src; zero_padding=zero_padding)
    return form_direct_image!(prepared)
end

"""Form a Shack--Hartmann detector-plane photon-rate mosaic explicitly."""
function form_shack_hartmann_rate!(wfs::ShackHartmannWFS, pupil,
    source=nothing)
    rate = shack_hartmann_rate_map(wfs, pupil, source)
    optics = prepare_wfs_optics(shack_hartmann_optics(wfs, source),
        pupil, rate)
    form_wfs_optical_products!(rate, pupil, optics)
    return rate
end

"""Form an SH photon-rate mosaic and acquire its detector observation."""
function acquire_shack_hartmann_observation!(wfs::ShackHartmannWFS,
    pupil, source, detector::Detector, rng)
    rate = form_shack_hartmann_rate!(wfs, pupil, source)
    frame = similar(rate.values)
    observation = WFSObservation(frame; units=:electron_count,
        layout=:lenslet_mosaic)
    acquisition = prepare_wfs_acquisition(detector, rate, observation;
        source=source)
    acquire_wfs_observation!(observation, rate, acquisition, rng)
    return observation
end

"""Form a Pyramid four-pupil detector-plane photon-rate product explicitly."""
function form_pyramid_rate!(wfs::PyramidWFS, pupil, source)
    front_end = PyramidOpticalFrontEnd(wfs, source)
    rate = pyramid_rate_map(front_end, pupil)
    optics = prepare_wfs_optics(front_end, pupil, rate)
    form_wfs_optical_products!(rate, pupil, optics)
    return rate
end

"""Form path-local Pyramid products for an asterism or extended source."""
function form_pyramid_path_rates!(wfs::PyramidWFS, pupils, source)
    front_end = PyramidOpticalFrontEnd(wfs, source)
    rates = pyramid_rate_map(front_end, pupils)
    optics = prepare_wfs_optics(front_end, pupils, rates)
    form_wfs_optical_products!(rates, pupils, optics)
    return rates
end

"""Form a Pyramid photon-rate frame and acquire its detector observation."""
function acquire_pyramid_observation!(wfs::PyramidWFS, pupil, source,
    detector::Detector, rng)
    rate = form_pyramid_rate!(wfs, pupil, source)
    observation = WFSObservation(similar(intensity_values(rate));
        units=:electron_count, layout=:four_pupil_mosaic)
    acquisition = prepare_wfs_acquisition(detector, rate, observation;
        source=source)
    acquire_wfs_observation!(observation, rate, acquisition, rng)
    return observation
end

"""Form a Bi-O-edge four-pupil detector-plane photon-rate product."""
function form_bi_o_edge_rate!(wfs::BiOEdgeWFS, pupil, source)
    front_end = BiOEdgeOpticalFrontEnd(wfs, source)
    rate = bi_o_edge_rate_map(front_end, pupil)
    optics = prepare_wfs_optics(front_end, pupil, rate)
    form_wfs_optical_products!(rate, pupil, optics)
    return rate
end

"""Form a Bi-O-edge photon-rate frame and acquire its detector observation."""
function acquire_bi_o_edge_observation!(wfs::BiOEdgeWFS, pupil, source,
    detector::Detector, rng)
    rate = form_bi_o_edge_rate!(wfs, pupil, source)
    observation = WFSObservation(similar(intensity_values(rate));
        units=:electron_count, layout=:four_pupil_mosaic)
    acquisition = prepare_wfs_acquisition(detector, rate, observation;
        source=source)
    acquire_wfs_observation!(observation, rate, acquisition, rng)
    return observation
end

function run_gpu_smoke_matrix(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    disable_scalar_backend!(B)
    failures = String[]
    rng = MersenneTwister(1)
    T = Float32
    atmosphere_step = T(1e-3)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    backend = AdaptiveOpticsSim.Backends.array_backend_selector(BackendArray)

    tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    lgs = LGSSource(; magnitude=0.0, wavelength=589e-9, altitude=90_000.0,
        laser_coordinates=(0.0, 0.0), photon_irradiance=one(T), T=T)
    spider_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
    apply_spiders!(spider_tel; thickness=0.5, angles_deg=[0.0, 90.0])
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
        Optics.intensity!(intensity, field)
        @assert intensity isa BackendArray

        amplitude = similar(prepared.wavefront.opd, T,
            tel.params.resolution, tel.params.resolution)
        fill!(amplitude, T(0.5))
        Optics.apply_amplitude!(field, amplitude, prepared.plan)
        Optics.intensity!(intensity, field)
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
        build_mask!(support, AnnularAperture(inner_radius=T(0.2), outer_radius=T(1), T=T); grid=Optics.default_mask_grid(support; T=T))
        apply_mask!(support, SpiderMask(thickness=T(0.1), angle_rad=T(pi / 4), T=T); grid=Optics.default_mask_grid(support; T=T))
        build_mask!(weighted, CircularAperture(radius=T(4), T=T); grid=Optics.pixel_mask_grid(weighted; T=T), inside=complex(T(inv(sqrt(2))), T(inv(sqrt(2)))))
        build_mask!(valid, SubapertureGridMask(threshold=T(0.1), T=T), support)
        @assert support isa BackendArray
        @assert weighted isa BackendArray
        @assert valid isa BackendArray
        @assert any(Array(valid))
        return valid
    end

    record_gpu_smoke!(failures, "detector_capture_none") do
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=2, T=T, backend=backend)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(acquisition; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "detector_capture_noise") do
        rate_map = gpu_direct_image(tel, src; zero_padding=2, T=T)
        det = Detector(noise=(NoisePhoton(), NoiseReadout(T(1e-3))), exposure_duration=1.0, qe=1.0,
            binning=2, background_flux=T(0.5), dark_current=T(0.1), T=T, backend=backend)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(acquisition; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "atmosphere_step") do
        atm = KolmogorovAtmosphere(tel; r0=0.2,
            reference_wavelength_m=T(500e-9), L0=25.0, T=T,
            backend=backend)
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
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        output = PupilFunction(tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, tel)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert atm.layers[1].generator.state.phase_rad isa BackendArray
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_multilayer_step_spiders") do
        atm = MultiLayerAtmosphere(spider_tel;
            r0=T(0.2),
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        output = PupilFunction(spider_tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, spider_tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(output, renderer, atm, epoch)
        @assert atm.layers[1].generator.state.phase_rad isa BackendArray
        @assert output.opd isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_multilayer_step") do
        atm = InfiniteMultiLayerAtmosphere(tel;
            r0=T(0.2),
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
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
        @assert atm.layers[1].screen.state.phase_rad isa BackendArray
        @assert atm.layers[1].screen.state.phase_scratch_rad isa BackendArray
        return output.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_multilayer_step_spiders") do
        atm = InfiniteMultiLayerAtmosphere(spider_tel;
            r0=T(0.2),
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
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
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
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
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
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
                reference_wavelength_m=T(500e-9),
                L0=T(25.0),
                fractional_cn2=T[0.7, 0.3],
                wind_speed=T[8.0, 4.0],
                wind_direction_deg=T[0.0, 90.0],
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
        atm = KolmogorovAtmosphere(tel; r0=0.2,
            reference_wavelength_m=T(500e-9), L0=25.0, T=T,
            backend=backend)
        ws = Atmospheres.PhaseStatsWorkspace(tel.params.resolution; T=T, backend=backend)
        screen, psd = Atmospheres.ft_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, return_psd=true)
        sh_screen = Atmospheres.ft_sh_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, subharmonics=true, n_levels=2, subharmonic_radius=1)
        rho = similar(atm.state.phase_rad, T, 4, 4)
        copyto!(rho, reshape(T[0.0, 0.02, 0.05, 0.1, 0.15, 0.25, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0], 4, 4))
        cov = Atmospheres.phase_covariance(rho, atm)
        freqs = similar(atm.state.freqs, T, 4)
        copyto!(freqs, T[0.1, 0.2, 0.3, 0.4])
        covmat = Atmospheres.covariance_matrix(freqs, freqs, atm)
        spectrum = Atmospheres.phase_spectrum(freqs, atm)
        freq_grid = similar(atm.state.phase_rad, T, 2, 2)
        copyto!(freq_grid, reshape(T[0.1, 0.2, 0.3, 0.4], 2, 2))
        spectrum_grid = Atmospheres.phase_spectrum(freq_grid, atm)
        @assert screen isa BackendArray
        @assert psd isa BackendArray
        @assert sh_screen isa BackendArray
        @assert cov isa BackendArray
        @assert covmat isa BackendArray
        @assert spectrum isa BackendArray
        @assert spectrum_grid isa BackendArray
        cov_ref = Atmospheres.phase_covariance(Array(rho), atm)
        covmat_ref = Atmospheres.covariance_matrix(Array(freqs), Array(freqs), atm)
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

    record_gpu_smoke!(failures, "shack_hartmann_photon_rate") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rate = form_shack_hartmann_rate!(wfs, pupil, src)
        @assert rate.values isa BackendArray
        return rate.values
    end

    record_gpu_smoke!(failures, "shack_hartmann_optics_and_acquisition") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=one(T), qe=one(T),
            binning=1, T=T, backend=backend)
        observation = acquire_shack_hartmann_observation!(wfs, pupil, src,
            det, rng)
        @assert observation_storage(observation) isa BackendArray
        return observation_storage(observation)
    end

    record_gpu_smoke!(failures, "shack_hartmann_spectral_rate_bundle") do
        bundle = SpectralBundle(fill(wavelength(src), 2), T[0.4, 0.6]; T=T)
        spectral = with_spectrum(src, bundle)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rates = form_shack_hartmann_rate!(wfs, pupil, spectral)
        @assert all(rate.values isa BackendArray for rate in rates)
        return rates[1].values
    end

    record_gpu_smoke!(failures, "shack_hartmann_spectral_distinct_wavelengths") do
        bundle = SpectralBundle(
            T[0.9 * wavelength(src), 1.1 * wavelength(src)],
            T[0.4, 0.6]; T=T)
        spectral = with_spectrum(src, bundle)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rates = form_shack_hartmann_rate!(wfs, pupil, spectral)
        @assert all(rate.values isa BackendArray for rate in rates)
        return rates[1].values
    end

    record_gpu_smoke!(failures, "shack_hartmann_extended_source_rate") do
        model = GaussianDiskSourceModel(sigma_arcsec=T(0.35), n_side=5, T=T)
        ext = with_extended_source(src, model)
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rate = form_shack_hartmann_rate!(wfs, pupil, ext)
        @assert rate.values isa BackendArray
        return rate.values
    end

    record_gpu_smoke!(failures, "shack_hartmann_spider_rate") do
        wfs = ShackHartmannWFS(spider_tel; n_lenslets=4, T=T, backend=backend)
        rate = form_shack_hartmann_rate!(wfs, spider_pupil, src)
        @assert rate.values isa BackendArray
        return rate.values
    end

    record_gpu_smoke!(failures, "shack_hartmann_lgs_rate") do
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rate = form_shack_hartmann_rate!(wfs, pupil, lgs)
        @assert rate.values isa BackendArray
        return rate.values
    end

    record_gpu_smoke!(failures, "shack_hartmann_asterism_rate") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        rate = form_shack_hartmann_rate!(wfs, pupil, ast)
        @assert rate.values isa BackendArray
        return rate.values
    end

    record_gpu_smoke!(failures, "shack_hartmann_asterism_acquisition") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmannWFS(tel; n_lenslets=4, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1, T=T, backend=backend)
        observation = acquire_shack_hartmann_observation!(wfs, pupil, ast,
            det, rng)
        @assert observation_storage(observation) isa BackendArray
        return observation_storage(observation)
    end

    record_gpu_smoke!(failures, "shack_hartmann_detector_equivalence") do
        cpu_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=CPUBackend())
        gpu_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        cpu_src = Source(band=:I, magnitude=0.0, T=T)
        gpu_src = Source(band=:I, magnitude=0.0, T=T)
        cpu_wfs = ShackHartmannWFS(cpu_tel; n_lenslets=4, T=T, backend=CPUBackend())
        gpu_wfs = ShackHartmannWFS(gpu_tel; n_lenslets=4, T=T, backend=backend)
        cpu_det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0,
            sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
        gpu_det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0,
            sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)

        cpu_pupil = PupilFunction(cpu_tel; T=T, backend=CPUBackend())
        gpu_pupil = PupilFunction(gpu_tel; T=T, backend=backend)
        cpu_observation = acquire_shack_hartmann_observation!(cpu_wfs,
            cpu_pupil, cpu_src, cpu_det, rng)
        gpu_observation = acquire_shack_hartmann_observation!(gpu_wfs,
            gpu_pupil, gpu_src, gpu_det, rng)

        cpu_export = Array(observation_storage(cpu_observation))
        gpu_export = Array(observation_storage(gpu_observation))

        @assert size(gpu_export) == size(cpu_export)
        @assert isapprox(gpu_export, cpu_export; rtol=1f-5, atol=1f-4)
        return gpu_export
    end

    record_gpu_smoke!(failures, "pyramid_photon_rate") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T, backend=backend)
        rate = form_pyramid_rate!(wfs, pupil, src)
        @assert intensity_values(rate) isa BackendArray
        return intensity_values(rate)
    end

    record_gpu_smoke!(failures, "pyramid_lgs_photon_rate") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T,
            backend=backend)
        rate = form_pyramid_rate!(wfs, pupil, lgs)
        @assert intensity_values(rate) isa BackendArray
        @assert all(isfinite, Array(intensity_values(rate)))
        return intensity_values(rate)
    end

    record_gpu_smoke!(failures, "pyramid_spectral_rate_bundle") do
        bundle = SpectralBundle(T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6]; T=T)
        poly = with_spectrum(src, bundle)
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T, backend=backend)
        rates = form_pyramid_rate!(wfs, pupil, poly)
        @assert all(rate -> intensity_values(rate) isa BackendArray, rates)
        return intensity_values(first(rates))
    end

    record_gpu_smoke!(failures, "pyramid_extended_source_rate_bundle") do
        model = GaussianDiskSourceModel(sigma_arcsec=T(0.35), n_side=5, T=T)
        ext = with_extended_source(src, model)
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T, backend=backend)
        pupils = ntuple(_ -> pupil, length(extended_source_asterism(ext)))
        rates = form_pyramid_path_rates!(wfs, pupils, ext)
        @assert all(rate -> intensity_values(rate) isa BackendArray, rates)
        return intensity_values(first(rates))
    end

    record_gpu_smoke!(failures, "pyramid_detector_acquisition") do
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1, T=T, backend=backend)
        observation = acquire_pyramid_observation!(wfs, pupil, src, det, rng)
        @assert observation_storage(observation) isa BackendArray
        return observation_storage(observation)
    end

    record_gpu_smoke!(failures, "bi_o_edge_photon_rate") do
        wfs = BiOEdgeWFS(tel; pupil_samples=4, modulation=0.0, T=T,
            backend=backend)
        rate = form_bi_o_edge_rate!(wfs, pupil, src)
        @assert intensity_values(rate) isa BackendArray
        @assert all(isfinite, Array(intensity_values(rate)))
        return intensity_values(rate)
    end

    record_gpu_smoke!(failures, "bi_o_edge_detector_acquisition") do
        wfs = BiOEdgeWFS(tel; pupil_samples=4, modulation=2.0, T=T,
            backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0,
            binning=1, T=T, backend=backend)
        observation = acquire_bi_o_edge_observation!(
            wfs, pupil, src, det, rng)
        @assert observation_storage(observation) isa BackendArray
        return observation_storage(observation)
    end

    record_gpu_smoke!(failures, "zernike_optics") do
        wfs = ZernikeWFS(tel; pupil_samples=4, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1, T=T, backend=backend)
        front_end = ZernikeOpticalFrontEnd(wfs, src)
        rate = zernike_rate_map(front_end, pupil)
        optics = prepare_wfs_optics(front_end, pupil, rate)
        observation = WFSObservation(similar(rate.values);
            units=:electron_count, layout=:zernike_pupil_image)
        acquisition = prepare_wfs_acquisition(det, rate, observation)
        form_wfs_optical_products!(rate, pupil, optics)
        acquire_wfs_observation!(observation, rate, acquisition, rng)
        @assert observation_storage(observation) isa BackendArray
        return observation_storage(observation)
    end

    record_gpu_smoke!(failures, "curvature_atmosphere_optics") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            reference_wavelength_m=T(500e-9),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction_deg=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=backend,
        )
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        atmospheric_pupil = PupilFunction(tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, tel, src)
        render_atmosphere!(atmospheric_pupil, renderer, atm, epoch)
        wfs = CurvatureWFS(tel; pupil_samples=4, T=T, backend=backend)
        front_end = CurvatureOpticalFrontEnd(wfs, src)
        rates = curvature_rate_maps(front_end, atmospheric_pupil)
        optics = prepare_wfs_optics(front_end, atmospheric_pupil, rates)
        form_wfs_optical_products!(rates, atmospheric_pupil, optics)
        @assert all(rate -> rate.values isa BackendArray, rates)
        return rates[1].values
    end

    record_gpu_smoke!(failures, "plant_step") do
        step_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        atm = KolmogorovAtmosphere(step_tel; r0=0.2,
            reference_wavelength_m=T(500e-9), L0=25.0, T=T,
            backend=backend)
        dm = DeformableMirror(step_tel; n_act=4, influence_width=0.3, T=T, backend=backend)
        wfs = ShackHartmannWFS(step_tel; n_lenslets=4, T=T, backend=backend)
        det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1, T=T, backend=backend)
        step_pupil = PupilFunction(step_tel; T=T, backend=backend)
        renderer = prepare_atmosphere_renderer(atm, step_tel, src)
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(step_pupil, renderer, atm, epoch)
        update_surface!(dm)
        apply_surface!(step_pupil, dm, DMAdditive())
        observation = acquire_shack_hartmann_observation!(wfs, step_pupil,
            src, det, rng)
        imaging = prepare_direct_imaging(step_pupil, src; zero_padding=2)
        rate_map = form_direct_image!(imaging)
        acquisition = prepare_detector_acquisition(det, rate_map)
        frame = capture!(acquisition; rng=rng)
        @assert observation_storage(observation) isa BackendArray
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "lift") do
        lift_tel = Telescope(resolution=16, diameter=8.0f0, central_obstruction=0.0f0, T=T, backend=backend)
        lift_src = Source(band=:I, magnitude=8.0, T=T)
        basis = backend_rand(B, T, 16, 16, 3)
        diversity = backend_zeros(B, T, 16, 16)
        model_opd = backend_zeros(B, T, 16, 16)
        forward = prepare_lift_forward_model(lift_tel, lift_src, basis,
            model_opd;
            diversity_opd=diversity, focal_resolution=32,
            zero_padding=2)
        rate_map = gpu_direct_image(lift_tel, lift_src;
            zero_padding=2, T=T)
        observation = LiFTObservation(forward, rate_map.values)
        specification = PR.LiFTSpecification(forward, observation)
        method = PR.LiFT(iterations=2,
            jacobian_method=PR.LiFTAnalyticJacobian(),
            solve_mode=PR.LiFTSolveNormalEquations(),
            mode_indices=(1, 2),
            model_scaling=PR.LiFTPhysicalRatePreservation(),
            check_convergence=false)
        execution = AdaptiveOpticsCalibration.KernelExecution(specification,
            KernelAbstractions.get_backend(rate_map.values); workgroup_size=64)
        plan = AdaptiveOpticsCalibration.prepare(method, execution)
        result = AdaptiveOpticsCalibration.allocate_result(plan)
        workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
        inputs = PR.LiFTInputs(observation.values)
        AdaptiveOpticsCalibration.process!(result, workspace, plan, inputs)
        coeffs = PR.lift_coefficients(result)
        @assert coeffs isa BackendArray
        return coeffs
    end

    if !isempty(failures)
        error("GPU smoke matrix failed:\n" * join(failures, "\n"))
    end

    println("gpu_smoke_matrix complete")
    return nothing
end

include("compiled_paths_helpers.jl")

@testset "Compiled fast filter batches" begin
    for T in (Float32, Float64), image_size in (5, 6), phantom_pad in (0, 2)
        fixture = compiled_paths_filtered_fixture(; image_size, phantom_pad)
        dm, opd, attenuation = (zeros(T, fixture.n, fixture.n) for _ in 1:3)
        prepare() = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(
            fixture.metadata, ref -> fixture.arrays[ref]; T,
            dm_opd=dm, entrance_opd=opd, entrance_attenuation=attenuation)
        owner, scalar = prepare(), prepare()
        fast = owner.workspace.fast
        @test size(fast.image.field) == (image_size, image_size, 8)
        @test size(fast.defocus.field) == (image_size+2phantom_pad, image_size+2phantom_pad, 4)
        @test collect(fast.defocus_slots) == [1, 2, 0, 0, 3, 4, 0, 0]
        retained = (fast.image.field, fast.image.transfer, fast.defocus.field, fast.defocus.transfer)
        rtol, atol = T === Float32 ? (3e-5, 3e-7) : (5e-12, 1e-14)
        for frame in 1:3
            for x in 1:fixture.n, y in 1:fixture.n
                dm[y,x] = T(frame*2e-8sin(0.2x+0.3y))
                opd[y,x] = T(frame*3e-8cos(0.4x-0.2y))
                attenuation[y,x] = T(frame*(0.01x+0.005y))
            end
            fill!(fast.image.field, complex(T(NaN)))
            fill!(fast.defocus.field, complex(T(NaN)))
            AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
            CP._enqueue_compiled_paths!(scalar, true, false)
            for name in (:fringed, :unfringed, :reference)
                @test getproperty(owner.products, name) ≈ getproperty(scalar.products, name) rtol=rtol atol=atol
            end
            for (actual, sources, free_ghost) in (
                    (owner.workspace.closed_envelopes, owner.workspace.closed_sources, true),
                    (owner.workspace.only_envelopes, owner.workspace.only_sources, false))
                expected = compiled_paths_fast_oracle(owner, sources, free_ghost)
                for (values, reference) in zip(actual, expected)
                    @test values ≈ reference rtol=rtol atol=atol
                end
            end
            @test all(((fast.image.field, fast.image.transfer, fast.defocus.field, fast.defocus.transfer) .=== retained))
        end
        @test compiled_paths_allocated(owner) == 0
    end
    # Execute the exact accelerator pack/unpack kernels on KA's CPU backend.
    # Odd/even grids and workgroups extending beyond their bounds exercise
    # every wrap direction, cropped access, padding and disabled ghost plane.
    ka = CP.KernelAbstractions
    for m in (5, 6), pad in (0, 2)
        p = m+2pad
        input = ComplexF64[complex(y, 2x-y) for y in 1:m, x in 1:m]
        field = fill(complex(NaN), p, p, 2)
        expected = zeros(ComplexF64, p, p)
        expected[pad+1:pad+m, pad+1:pad+m] .= input
        for enabled in (true, false)
            CP._compiled_paths_pack_filter_kernel!(ka.CPU(), (4, 4))(
                field, input, 2, pad, enabled; ndrange=(p, p))
            ka.synchronize(ka.CPU())
            @test field[:,:,2] == (enabled ? FFTW.ifftshift(expected) : zeros(ComplexF64,p,p))
            @test all(isnan, field[:,:,1])
        end
        defocus = reshape(ComplexF64.(1:2p^2), p, p, 2)
        image = zeros(ComplexF64, m, m, 2)
        CP._compiled_paths_pack_defocus_kernel!(ka.CPU(), (4, 4))(
            image, defocus, 1, 2, pad; ndrange=(m, m))
        ka.synchronize(ka.CPU())
        @test image[:,:,1] == FFTW.ifftshift(FFTW.fftshift(defocus[:,:,2])[pad+1:pad+m,pad+1:pad+m])
        phase = cis.(0.03 .* real.(input))
        out = similar(input)
        CP._compiled_paths_unpack_filter_kernel!(ka.CPU(), (4, 4))(
            out, image, phase, 1; ndrange=(m, m))
        ka.synchronize(ka.CPU())
        @test out == FFTW.fftshift(image[:,:,1]) .* phase
    end
end

@testset "Frozen compiled_paths native execution" begin
    for forward in (true,false), from_waist in (true,false)
        fixture = compiled_paths_synthetic_fixture(;forward,from_waist)
        dm = [1e-8*sin(0.7x+0.3y) for y in 1:fixture.n,x in 1:fixture.n]
        entrance_opd = [2e-8*cos(0.2x-0.4y) for y in 1:fixture.n,x in 1:fixture.n]
        attenuation = [0.01x+0.005y for y in 1:fixture.n,x in 1:fixture.n]
        owner = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata, ref->fixture.arrays[ref];
            dm_opd=dm, entrance_opd, entrance_attenuation=attenuation)
        @test owner.inputs.dm_opd === dm
        @test owner.inputs.entrance_opd === entrance_opd
        @test owner.inputs.entrance_attenuation === attenuation
        @test all(iszero,owner.products.fringed)
        products = AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
        expected = compiled_paths_synthetic_oracle(fixture,dm,entrance_opd,attenuation)
        for name in (:fringed,:unfringed,:reference)
            @test getproperty(products,name) ≈ getproperty(expected,name) rtol=2e-12 atol=1e-15
            @test all(isfinite,getproperty(products,name))
        end
        before = copy(products.fringed)
        dm[3,5] += 2e-7
        AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
        expected = compiled_paths_synthetic_oracle(fixture,dm,entrance_opd,attenuation)
        @test products.fringed ≈ expected.fringed rtol=2e-12 atol=1e-15
        @test norm(products.fringed-before) > 1e-5
        @test compiled_paths_allocated(owner) == 0
    end
    fixture = compiled_paths_synthetic_fixture()
    caller_owned = CP._prepare_compiled_paths(fixture.metadata,
        ref->fixture.arrays[ref], nothing)
    @test caller_owned.execution_context === nothing
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(caller_owned)
    @test CP._enqueue_compiled_paths!(caller_owned) === caller_owned.products
    standalone = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata, ref->fixture.arrays[ref])
    AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(standalone)
    @test caller_owned.products.fringed ≈ standalone.products.fringed
    @test compiled_paths_allocated(standalone) == 0
    dark=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,ref->fixture.arrays[ref];
        entrance_attenuation=fill(1000.0,fixture.n,fixture.n))
    @test all(iszero,AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(dark).fringed)
    # Execute the exact device reduction kernels on KA's CPU backend. This
    # checks fixed-block padding, deterministic tree order and dynamic power
    # without requiring an accelerator in ordinary package tests.
    ka = CP.KernelAbstractions
    # Run the accelerator copy kernels on irregular CPU grids. KA masks the
    # extra lanes in the final workgroup; the prepared shape and padding
    # contracts prove each explicitly inbounds source access.
    for shape in ((1,1), (3,5), (4,6), (17,19))
        input = reshape(Float64.(1:prod(shape)), shape)
        out = similar(input)
        for inverse in (false,true)
            sy = inverse ? fld(shape[1],2) : cld(shape[1],2)
            sx = inverse ? fld(shape[2],2) : cld(shape[2],2)
            CP._compiled_paths_shift_kernel!(ka.CPU(),(4,4))(
                out,input,sy,sx; ndrange=shape)
            ka.synchronize(ka.CPU())
            @test out == (inverse ? FFTW.ifftshift(input) : FFTW.fftshift(input))
        end
        crop_shape = max.(1,shape .- 1)
        cropped = zeros(crop_shape)
        for start in (0,minimum(shape .- crop_shape))
            CP._compiled_paths_crop_kernel!(ka.CPU(),(4,4))(
                cropped,input,start; ndrange=crop_shape)
            ka.synchronize(ka.CPU())
            @test cropped == input[1+start:crop_shape[1]+start,1+start:crop_shape[2]+start]
        end
        for pad in (0,2)
            padded = zeros(shape .+ 2pad)
            CP._compiled_paths_pad_kernel!(ka.CPU(),(4,4))(
                padded,input,pad; ndrange=size(padded))
            ka.synchronize(ka.CPU())
            expected = zeros(size(padded))
            expected[1+pad:shape[1]+pad,1+pad:shape[2]+pad] .= input
            @test padded == expected
        end
    end
    # The largest case has more than 256 partials, exercising each lane's
    # repeated addition in the second reduction stage.
    for shape in ((3,3), (16,17), (33,35), (257,259))
        field = ComplexF64[cis(0.17x+0.21y)*(1+0.002x*y) for y in 1:shape[1], x in 1:shape[2]]
        partials = zeros(Float64,cld(length(field),256))
        power = zeros(Float64,1)
        expected = sum(abs2,field)
        first_run = copy(field)
        for _ in 1:2
            copyto!(field,first_run)
            CP._compiled_paths_power_partials_kernel!(ka.CPU(),256)(
                partials,field,length(field); ndrange=256length(partials))
            CP._compiled_paths_power_finalize_kernel!(ka.CPU(),256)(
                power,partials,length(partials); ndrange=256)
            ka.synchronize(ka.CPU())
            @test power[1] ≈ expected rtol=8eps(Float64)
            CP._compiled_paths_normalize_kernel!(ka.CPU(),256)(
                field,power,length(field); ndrange=length(field))
            ka.synchronize(ka.CPU())
            @test sum(abs2,field) ≈ 1 rtol=8eps(Float64)
        end
    end
    # The compiler calibration is centred; the shared AOS detector filter
    # consumes corner-order frequencies. Use a nonseparable anisotropic MTF
    # on an odd padded grid to expose shift or axis mistakes.
    mtf_metadata = deepcopy(fixture.metadata)
    mtf_metadata["fast"]["apply_pixel_mtf"] = true
    mtf_arrays = deepcopy(fixture.arrays)
    f = FFTW.fftshift(FFTW.fftfreq(7))
    mtf_arrays["pixel"] = ComplexF64[
        exp(-3fy^2-2fx^2)*(1+0.15sin(2pi*fy)*sin(2pi*fx)) for fy in f,fx in f]
    filtered = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(mtf_metadata, ref->mtf_arrays[ref])
    AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(filtered)
    @test standalone.workspace.pixel_mtf === nothing
    @test !hasproperty(filtered.workspace,:pixel_fft)
    for (mtf_owner,name) in zip(filtered.workspace.pixel_mtf,(:fringed,:unfringed,:reference))
        @test mtf_owner.input === getproperty(filtered.products,name)
        @test !Base.mightalias(mtf_owner.input,mtf_owner.products.output)
        @test mtf_owner.plan.transfer == real.(FFTW.ifftshift(mtf_arrays["pixel"]))
        raw = getproperty(standalone.products,name)
        padded = zeros(7,7)
        padded[2:6,2:6] .= raw
        expected = real.(FFTW.ifft(FFTW.fft(padded) .* real.(FFTW.ifftshift(mtf_arrays["pixel"]))))[2:6,2:6]
        @test getproperty(filtered.products,name) ≈ expected rtol=2e-12 atol=1e-15
    end
    @test compiled_paths_allocated(filtered) == 0
    mtf_arrays["pixel"][1,1] += 0.01im
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(mtf_metadata, ref->mtf_arrays[ref])
    for name in ("phantom_pad","pixel_pad")
        bad = deepcopy(fixture.metadata)
        bad["fast"][name] = -1
        @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->fixture.arrays[ref])
        max_pad = (typemax(Int) - fixture.m) ÷ 2
        for oversized in (max_pad + 1, typemax(Int))
            bad["fast"][name] = oversized
            @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->fixture.arrays[ref])
        end
    end
    # typemax padding previously wrapped the padded size to 3, permitting a
    # matching defocus transfer and an unsafe crop from an undersized buffer.
    bad = deepcopy(fixture.metadata)
    bad_arrays = deepcopy(fixture.arrays)
    bad["fast"]["phantom_pad"] = typemax(Int)
    bad["fast_paths"][1]["defocus_transfer"] = ["wrapped_defocus"]
    bad_arrays["wrapped_defocus"] = ones(ComplexF64,3,3)
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->bad_arrays[ref])
    bad = deepcopy(fixture.metadata)
    bad["schema_version"] = 99
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->fixture.arrays[ref])
    bad = deepcopy(fixture.metadata)
    bad["array_orientation"] = "xy"
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->fixture.arrays[ref])
    bad = deepcopy(fixture.metadata)
    bad["propagations"][1]["prefix"][1]["input"] = "unknown"
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad, ref->fixture.arrays[ref])
    @test_throws DimensionMismatch AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
        ref->fixture.arrays[ref]; dm_opd=zeros(3,3))
    for value in (NaN, Inf)
        bad_arrays = deepcopy(fixture.arrays)
        bad_arrays["entrance"][1,1] = value
        @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata, ref->bad_arrays[ref])
        for role in ("basis", "coefficient", "visibility")
            bad_metadata = deepcopy(fixture.metadata)
            bad_arrays = deepcopy(fixture.arrays)
            bad_arrays["nonfinite"] = fill(ComplexF64(value), 1, 1)
            bad_metadata["fast"][role] = role == "visibility" ? ["nonfinite"] : "nonfinite"
            @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(bad_metadata, ref->bad_arrays[ref])
        end
    end
    overflow_arrays = deepcopy(fixture.arrays)
    overflow_arrays["entrance"][1,1] = 1e100
    @test_throws ArgumentError AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
        ref->overflow_arrays[ref]; T=Float32)
end

@testset "Compiled optical path artifact and live grids" begin
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory)
        f=model.fixture
        dm=[1e-8*sin(0.7x+0.3y) for y in 1:8,x in 1:8]
        opd=[2e-8*cos(0.2x-0.4y) for y in 1:8,x in 1:8]
        attenuation=[0.01x+0.005y for y in 1:8,x in 1:8]
        graph=compiled_paths_test_graph(model;dm,opd,amplitude=attenuation)
        step_graph!(graph)
        expected=compiled_paths_synthetic_oracle(f,dm,opd,attenuation)
        @test graph_output(graph,:fringed)≈expected.fringed rtol=2e-12 atol=1e-15
        @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(graph,:optics))==0.001
        # Different physical pitch and support despite a familiar square shape.
        coarse=(size=(4,4),axes=(1,2),sampling_m=(2f.dx,2f.dx),center_pixel=(3.0,3.0))
        dm4=[1e-7*(2y-x) for y in 1:4,x in 1:4]
        coarse_graph=compiled_paths_test_graph(model;dm=dm4,opd=zeros(4,4),amplitude=zeros(4,4),grid=coarse)
        step_graph!(coarse_graph)
        native=AlgorithmGraphs.prepared_graph_node(coarse_graph,:optics).propagation.inputs.dm_opd
        expected_native=[sum(max(0,1-abs((y-5)/2+3-i))*max(0,1-abs((x-5)/2+3-j))*dm4[i,j]
            for i in 1:4,j in 1:4) for y in 1:8,x in 1:8]
        @test native≈expected_native
        native_graph=compiled_paths_test_graph(model;dm=expected_native)
        step_graph!(native_graph)
        @test graph_output(coarse_graph,:fringed)≈graph_output(native_graph,:fringed)
        # A reference-zero sample owns baked opacity; additional gain there is rejected.
        reference=ones(8,8);reference[2,3]=0
        amplitude=copy(reference)
        ratio_graph=compiled_paths_test_graph(model;mode="reference_ratio",reference,dm,opd,amplitude)
        step_graph!(ratio_graph)
        @test graph_output(ratio_graph,:fringed)≈compiled_paths_synthetic_oracle(f,dm,opd,zeros(8,8)).fringed
        ratio_owner=AlgorithmGraphs.prepared_graph_node(ratio_graph,:optics)
        @test all(iszero,ratio_owner.propagation.inputs.entrance_attenuation)
        @test !hasproperty(ratio_owner.propagation.inputs,:entrance_opacity)
        # Changed support is rejected rather than encoded through an unqualified
        # product of independently interpolated opacity and attenuation.
        amplitude[3,4]=0
        @test_throws ArgumentError step_graph!(ratio_graph)
        amplitude[3,4]=reference[3,4]
        reset_graph!(ratio_graph);step_graph!(ratio_graph)
        @test all(isfinite,graph_output(ratio_graph,:fringed))
        amplitude[2,3]=0.5
        @test_throws ArgumentError step_graph!(ratio_graph)
        amplitude.=reference
        reset_graph!(ratio_graph);step_graph!(ratio_graph)
        @test graph_output(ratio_graph,:fringed)≈compiled_paths_synthetic_oracle(f,dm,opd,zeros(8,8)).fringed
        # Cold I/O rejects altered hashes, shape contracts and escaping assets.
        @test_throws ArgumentError AlgorithmGraphs._read_compiled_optical_model(model.path,"0"^64)
        metadata=deepcopy(model.metadata)
        metadata["arrays"]["entrance"]["path"]="../outside"
        path=joinpath(directory,"escape.yaml");write(path,YAML.write(metadata))
        _,reader=AlgorithmGraphs._read_compiled_optical_model(path,bytes2hex(sha256(read(path))))
        @test_throws ArgumentError reader("entrance")
        FITS(joinpath(directory,"entrance"),"r+") do file
            write(file[1],ones(8,8))
        end
        @test_throws ArgumentError compiled_paths_test_graph(model)
    end
end

@testset "Compiled output sampling calibration preserves samples" begin
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory)
        for T in (Float32,Float64)
            dm=T[1e-8*sin(x+2y) for y in 1:8,x in 1:8]
            opd=zeros(T,8,8)
            amplitude=zeros(T,8,8)
            native=compiled_paths_test_graph(model;dm,opd,amplitude)
            calibration=(sampling_m=15e-6,provenance="synthetic detector coordinate calibration")
            calibrated=compiled_paths_test_graph(model;dm,opd,amplitude,
                output_sampling_calibration=calibration)
            owner=AlgorithmGraphs.prepared_graph_node(calibrated,:optics)
            @test optical_sampling_m(owner)===T(calibration.sampling_m)
            @test owner.output_sampling_calibration.provenance==calibration.provenance
            @test AlgorithmGraphs.prepared_graph_node(native,:optics).output_sampling_calibration===nothing
            step_graph!(native);step_graph!(calibrated)
            for name in (:fringed,:unfringed,:reference)
                @test graph_output(calibrated,name)==graph_output(native,name)
                @test size(graph_output(calibrated,name))==(5,5)
            end
        end
        for calibration in ((sampling_m=0,provenance="source"),
                (sampling_m=-1e-6,provenance="source"),
                (sampling_m=Inf,provenance="source"),
                (sampling_m=NaN,provenance="source"),
                (sampling_m=true,provenance="source"),
                (sampling_m=(15e-6,15e-6),provenance="source"),
                (sampling_m=15e-6,provenance=" "),
                (sampling_m=15e-6,provenance=17),
                (sampling_m=15e-6,), (provenance="source",),
                (sampling_m=15e-6,provenance="source",roi_origin=(1,1)), 15e-6)
            @test_throws AlgorithmGraphError compiled_paths_test_graph(model;
                output_sampling_calibration=calibration)
        end
        for pitch in (1e-300,1e300)
            @test_throws AlgorithmGraphError compiled_paths_test_graph(model;
                dm=zeros(Float32,8,8),opd=zeros(Float32,8,8),amplitude=zeros(Float32,8,8),
                output_sampling_calibration=(sampling_m=pitch,provenance="source"))
        end

        # The ordinary observation validator consumes the calibrated coordinates;
        # matching dimensions alone cannot admit a mismatched detector grid.
        detector=ingaas_detector_acquisition_node(:detector;rows=5,columns=5,
            pixel_sampling_m=(15e-6,15e-6),spectral_channel=:synthetic,
            exposure_duration_s=1e-3,quantum_efficiency=0.5,photon_noise=false,
            readout_noise=false,rng_seed=17,photon_rate_schema="test.rate/1",
            frame_schema="test.frame/1",T=Float64)
        response=prepare_algorithm_graph(algorithm_graph((detector,);
            inputs=(graph_input(:photon_rate,:detector=>:photon_rate,zeros(5,5)),),
            outputs=(graph_output(:frame,:detector=>:frame),)))
        specification=(graph=response,input=:photon_rate,output=:frame,
            normalized_to_photon_rate=1e6,pixel_sampling_m=(15e-6,15e-6),
            exposure_duration_s=1e-3,frame_schema="test.frame/1")
        native=compiled_paths_test_graph(model)
        observations=(sensor=(output=:fringed,sampling_node=:optics,response=specification),)
        @test_throws AlgorithmGraphError prepare_external_pupil((sensor=native,);
            exposure_durations=(sensor=ModelDuration(1_000_000),),observations,
            input_names=(opd=:opd,amplitude=:amplitude))
        calibrated=compiled_paths_test_graph(model;output_sampling_calibration=(
            sampling_m=15e-6,provenance="synthetic detector coordinate calibration"))
        prepared=prepare_external_pupil((sensor=calibrated,);
            exposure_durations=(sensor=ModelDuration(1_000_000),),observations,
            input_names=(opd=:opd,amplitude=:amplitude))
        images=Matrix{Float64}[]
        emit(name,image,metadata)=push!(images,copy(image))
        @test step_external_pupil!(emit,prepared,zeros(8,8),zeros(8,8),
            ModelTimestamp(0),ModelDuration(1_000_000);domain=ntuple(_->UInt8(1),16),
            generation=UInt64(1),sequence=UInt64(1))
        @test only(images)≈500 .* graph_output(calibrated,:fringed)

        # Shared scalar Proper prefixes must never replace compiled/fitted
        # propagation, even when multiple calibrated paths are observed.
        other=compiled_paths_test_graph(model;output_sampling_calibration=(
            sampling_m=15e-6,provenance="synthetic detector coordinate calibration"))
        compiled=prepare_external_pupil((;calibrated,other);
            exposure_durations=(calibrated=ModelDuration(1_000_000),other=ModelDuration(1_000_000)),
            observations=(calibrated=(output=:fringed,sampling_node=:optics),
                other=(output=:reference,sampling_node=:optics)),
            input_names=(opd=:opd,amplitude=:amplitude))
        @test compiled.propagation===nothing
        @test AlgorithmGraphs.optical_execution_graph(compiled,:calibrated)===calibrated
        @test AlgorithmGraphs.optical_execution_graph(compiled,:other)===other
        @test step_external_pupil!((args...)->nothing,compiled,zeros(8,8),zeros(8,8),
            ModelTimestamp(0),ModelDuration(1_000_000);domain=ntuple(_->UInt8(1),16),
            generation=UInt64(1),sequence=UInt64(1))
        @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(calibrated,:optics))==15e-6
    end
end

@testset "Compiled source-mode live input contract" begin
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory;source_mode=true)
        dm=zeros(8,8);opd=zeros(8,8);amplitude=ones(8,8);reference=ones(8,8)
        graph=compiled_paths_test_graph(model;mode="reference_ratio",reference,dm,opd,amplitude)
        step_graph!(graph)
        baseline=copy(graph_output(graph,:fringed))
        opd.=[1e-7sin(i+2j) for i in 1:8,j in 1:8]
        amplitude.=0.8
        step_graph!(graph)
        @test graph_output(graph,:fringed)==baseline # Frozen internal-source entrance.
        for input in (opd,amplitude)
            saved=input[1,1]
            input[1,1]=NaN
            @test_throws ArgumentError step_graph!(graph)
            input[1,1]=saved
            reset_graph!(graph)
        end
        amplitude[1,1]=-0.1
        @test_throws ArgumentError step_graph!(graph)
        amplitude[1,1]=0.8
        reset_graph!(graph);step_graph!(graph)
        @test graph_output(graph,:fringed)==baseline
    end
end

@testset "Compiled input validation and ratio device kernels on KA CPU" begin
    ka=CP.KernelAbstractions
    for T in (Float32,Float64)
        reference=ones(T,3,5);reference[1,1]=0
        amplitude=T[.2+.03i+.02j for i in 1:3,j in 1:5];amplitude[1,1]=0
        expected=similar(amplitude);actual=similar(amplitude)
        failure=zeros(UInt32,6)
        CP._compiled_paths_ratio!(expected,amplitude,reference)
        CP._compiled_paths_validate_amplitude_kernel!(ka.CPU(),32)(failure,amplitude,reference;ndrange=length(amplitude))
        CP._compiled_paths_ratio_kernel!(ka.CPU(),32)(actual,amplitude,reference,failure;ndrange=length(actual))
        ka.synchronize(ka.CPU())
        @test actual≈expected rtol=4eps(T) atol=eps(T)
        @test all(iszero,failure)
        for (a,r,slot) in ((T(NaN),one(T),3),(T(-.1),one(T),3),
                (T(.2),zero(T),4),(zero(T),one(T),5))
            fill!(failure,0);amplitude[2,3]=a;reference[2,3]=r
            CP._compiled_paths_validate_amplitude_kernel!(ka.CPU(),32)(failure,amplitude,reference;ndrange=length(amplitude))
            ka.synchronize(ka.CPU())
            @test findall(!iszero,failure)==[slot]
            @test_throws ArgumentError CP._compiled_paths_validate_amplitude!(amplitude,reference)
        end
        amplitude[2,3]=floatmax(T);reference[2,3]=floatmin(T);fill!(failure,0)
        CP._compiled_paths_ratio_kernel!(ka.CPU(),32)(actual,amplitude,reference,failure;ndrange=length(actual))
        ka.synchronize(ka.CPU())
        @test findall(!iszero,failure)==[6]
        @test all(isfinite,actual)
        @test_throws ArgumentError CP._compiled_paths_ratio!(expected,amplitude,reference)
        for slot in (1,2), value in (T(NaN),T(Inf))
            fill!(failure,0);input=zeros(T,3,5);input[end]=value
            CP._compiled_paths_validate_finite_kernel!(ka.CPU(),32)(failure,input,slot;ndrange=length(input))
            ka.synchronize(ka.CPU())
            @test findall(!iszero,failure)==[slot]
        end
    end
end

@testset "Compiled optical common prefixes" begin
    for summation in (:path_pairs, :coherent_modes)
        fixture = compiled_paths_shared_fixture()
        dm = zeros(fixture.n, fixture.n)
        opd, attenuation = copy(dm), copy(dm)
        owner = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(
            fixture.metadata, ref -> fixture.arrays[ref];
            dm_opd=dm, entrance_opd=opd, entrance_attenuation=attenuation, summation)
        @test length(owner.workspace.shared_prefixes) == 2
        @test all(p -> p.count == 5, owner.plan.shared_prefixes)
        @test [p.cache_index for p in owner.plan.shared_prefixes] == [1, 2, 1, 2, 1, 2]
        oracle = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(
            fixture.metadata, ref -> fixture.arrays[ref];
            dm_opd=dm, entrance_opd=opd, entrance_attenuation=attenuation, summation)
        before = nothing
        for frame in 1:5
            @. dm = 8e-8 * sin(frame * 0.3 * (1:fixture.n))
            @. opd = 2e-8 * cos(frame * 0.7 * (1:fixture.n))
            @. attenuation = frame * 0.02 * (1:fixture.n)
            # Poison retained fields to prove every invocation refreshes them.
            foreach(field -> fill!(field, complex(NaN)), owner.workspace.shared_prefixes)
            CP._enqueue_compiled_paths!(oracle, false)
            AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
            for name in (:fringed, :unfringed, :reference)
                @test getproperty(owner.products, name) == getproperty(oracle.products, name)
            end
            if before !== nothing
                @test owner.products.fringed != before
            end
            before = copy(owner.products.fringed)
        end
        @test compiled_paths_allocated(owner) == 0
    end
    fixture = compiled_paths_shared_fixture()
    fixture.arrays["other_entrance"] = copy(fixture.arrays["entrance"])
    # Different entrance bindings and different leading live wavelengths must
    # not share, even if another spectral anchor still has a common prefix.
    fixture.metadata["propagations"][3]["entrance"] = "other_entrance"
    fixture.metadata["propagations"][5]["prefix"][1]["lam"] *= 1.01
    owner = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(
        fixture.metadata, ref -> fixture.arrays[ref])
    @test [p.count for p in owner.plan.shared_prefixes] == [0, 5, 0, 5, 0, 5]
    independent = AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(
        fixture.metadata, ref -> fixture.arrays[ref])
    CP._enqueue_compiled_paths!(independent, false)
    AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
    @test owner.products.fringed == independent.products.fringed
    mktempdir() do directory
        model = compiled_paths_yaml_fixture(directory; fixture=compiled_paths_shared_fixture())
        graph = compiled_paths_test_graph(model)
        step_graph!(graph)
        original = copy(graph_output(graph, :fringed))
        owner = AlgorithmGraphs.prepared_graph_node(graph, :optics)
        owner.inputs.dm_opd .= 1e-7
        step_graph!(graph)
        reset_graph!(graph)
        @test all(iszero, graph_output(graph, :fringed))
        fill!(owner.inputs.dm_opd, 0)
        foreach(field -> fill!(field, complex(NaN)), owner.propagation.workspace.shared_prefixes)
        step_graph!(graph)
        @test graph_output(graph, :fringed) == original
    end
end

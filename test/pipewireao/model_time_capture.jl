using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.AlgorithmGraphs: ModelDuration, ModelTimestamp
using PipeWireAO
using Test

function pipewire_acquisition_fixture(storage)
    native = PipeWireAO.LibPipeWire
    metas = native.spa_meta[
        native.spa_meta(
            native.SPA_META_Acquisition,
            UInt32(sizeof(native.spa_meta_acquisition)),
            Ptr{Cvoid}(Base.pointer(storage)),
        ),
    ]
    spa_buffer = Ref(native.spa_buffer(
        UInt32(1),
        UInt32(0),
        pointer(metas),
        C_NULL,
    ))
    native_buffer = Ref(native.pw_buffer(
        Base.unsafe_convert(Ptr{native.spa_buffer}, spa_buffer),
        C_NULL,
        UInt64(0),
        UInt64(0),
        UInt64(0),
    ))
    buffer = StreamBuffer(
        Base.unsafe_convert(Ptr{native.pw_buffer}, native_buffer),
    )
    return (; buffer, metas, spa_buffer, native_buffer)
end

@testset "PipeWireAO acquisition model-time capture" begin
    native = PipeWireAO.LibPipeWire
    words = div(Int(native.SPA_META_ACQUISITION_SIZE), sizeof(UInt64))
    origin_storage = zeros(UInt64, words)
    sample_storage = zeros(UInt64, words)
    origin_fixture = pipewire_acquisition_fixture(origin_storage)
    sample_fixture = pipewire_acquisition_fixture(sample_storage)

    domain = AcquisitionDomain(ntuple(
        index -> index == 1 ? UInt8(1) : UInt8(0),
        Val(ACQUISITION_DOMAIN_SIZE),
    ))
    reference = AcquisitionPtpReference(
        PtpClockIdentity(ntuple(
            index -> index == 1 ? UInt8(2) : UInt8(0),
            Val(ACQUISITION_PTP_CLOCK_ID_SIZE),
        )),
        7,
    )

    GC.@preserve origin_storage sample_storage origin_fixture sample_fixture begin
        origin_metadata = buffer_acquisition(origin_fixture.buffer)
        sample_metadata = buffer_acquisition(sample_fixture.buffer)
        initialize_acquisition!(origin_metadata)
        initialize_acquisition!(sample_metadata)
        set_acquisition_identity!(
            origin_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(81)),
        )
        set_acquisition_identity!(
            sample_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(82)),
        )
        set_acquisition_exposure_start_ptp!(
            origin_metadata,
            1_000_000,
            10,
            reference,
        )
        set_acquisition_exposure_start_ptp!(
            sample_metadata,
            1_001_020,
            12,
            reference,
        )
        set_acquisition_exposure_duration!(origin_metadata, 1_000)
        set_acquisition_exposure_duration!(sample_metadata, 1_000)

        origin = capture_model_time_origin(origin_metadata)
        first_capture = capture_model_timestamp(origin_metadata, origin)
        second_capture = capture_model_timestamp(sample_metadata, origin)
        @test model_timestamp(first_capture) == ModelTimestamp(0)
        @test model_timestamp(second_capture) == ModelTimestamp(1_020)
        @test model_time_uncertainty(second_capture) == ModelDuration(12)
        provenance = model_time_provenance(second_capture)
        @test provenance.identity.sequence == UInt64(82)
        @test provenance.exposure_start_nanoseconds == Int64(1_001_020)
        @test provenance.timebase == ACQUISITION_TIMEBASE_TAI
        @test provenance.ptp_reference == reference
        @test provenance.exposure_duration_nanoseconds == UInt64(1_000)
        @test isbitstype(typeof(provenance))

        replay = prepare_captured_model_time_driver((
            first_capture,
            second_capture,
        ))
        @test advance_model_time!(replay) == ModelTimestamp(0)
        @test advance_model_time!(replay) == ModelTimestamp(1_020)
        @test model_time_exhausted(replay)

        capture_model_timestamp(sample_metadata, origin)
        @test @allocated(capture_model_timestamp(
            sample_metadata,
            origin,
        )) == 0

        initialize_acquisition!(sample_metadata)
        set_acquisition_identity!(
            sample_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(83)),
        )
        set_acquisition_exposure_start_ptp!(
            sample_metadata,
            999_999,
            12,
            reference,
        )
        set_acquisition_exposure_duration!(sample_metadata, 1_000)
        @test_throws AlgorithmGraphError capture_model_timestamp(
            sample_metadata,
            origin,
        )

        initialize_acquisition!(sample_metadata)
        set_acquisition_identity!(
            sample_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(84)),
        )
        set_acquisition_exposure_start!(sample_metadata, 1_002_000, 12)
        set_acquisition_exposure_duration!(sample_metadata, 1_000)
        @test_throws AlgorithmGraphError capture_model_timestamp(
            sample_metadata,
            origin,
        )

        other_reference = AcquisitionPtpReference(
            PtpClockIdentity(ntuple(
                index -> index == 1 ? UInt8(3) : UInt8(0),
                Val(ACQUISITION_PTP_CLOCK_ID_SIZE),
            )),
            7,
        )
        initialize_acquisition!(sample_metadata)
        set_acquisition_identity!(
            sample_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(85)),
        )
        set_acquisition_exposure_start_ptp!(
            sample_metadata,
            1_002_000,
            12,
            other_reference,
        )
        set_acquisition_exposure_duration!(sample_metadata, 1_000)
        @test_throws AlgorithmGraphError capture_model_timestamp(
            sample_metadata,
            origin,
        )

        initialize_acquisition!(sample_metadata)
        set_acquisition_identity!(
            sample_metadata,
            AcquisitionIdentity(domain, UInt64(4), UInt64(86)),
        )
        set_acquisition_exposure_start_ptp!(
            sample_metadata,
            1_002_000,
            12,
            reference,
        )
        @test_throws AlgorithmGraphError capture_model_timestamp(
            sample_metadata,
            origin,
        )
    end
end

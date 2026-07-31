@testset "CCD detector" begin
    zero_psf = zeros(4, 4)
    det_ccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    frame_ccd_cic = capture!(det_ccd_cic, zero_psf; rng=MersenneTwister(11))
    @test sum(frame_ccd_cic) > 0
    @test supports_clock_induced_charge(det_ccd_cic.params.sensor)
    det_ccd_cic_long = Detector(integration_time=10.0, noise=NoiseNone(),
        qe=1.0, sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    @test capture!(det_ccd_cic_long, zero_psf; rng=MersenneTwister(11)) ==
        frame_ccd_cic
    @test_throws InvalidConfiguration CCDSensor(
        clock_induced_charge_per_frame=-1.0)

    det_default_ccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, sensor=CCDSensor())
    @test detector_export_metadata(det_default_ccd).frame_response == :none
end

# API Reference

This is a guide to the exported public API. It is organized by subsystem rather
than by source file.

## Core types and utilities

- `Workspace`, `ensure_psf_buffers!`
- `FidelityProfile`, `ScientificProfile`, `FastProfile`,
  `default_fidelity_profile`
- `Telemetry`, `TelemetryRow`, `record!`
- `ClosedLoopTrace`, `ClosedLoopTraceRow`
- `GSCClosedLoopTrace`, `GSCClosedLoopTraceRow`
- `GSCAtmosphereReplayTrace`, `GSCAtmosphereReplayTraceRow`
- `snapshot_config`, `write_config_toml`, `write_config_json`
- `write_telemetry_csv`
- `bin2d`, `poisson_noise!`, `poisson_sample`
- `ParallelConfig`, `with_parallel_config`
- `InversePolicy`, `ExactPseudoInverse`, `TSVDInverse`, `TikhonovInverse`
- Exceptions: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`

## Optical elements

- `Telescope`, `TelescopeParams`, `TelescopeState`
- `generate_pupil!`, `reset_opd!`, `apply_opd!`, `set_pupil!`, `apply_spiders!`
- `Source`, `SourceParams`, `LGSSource`, `LGSSourceParams`, `wavelength`
- `Asterism`, `coordinates_xy_arcsec`, `compute_psf!`,
  `psf_pixel_scale_arcsec`, `ensure_psf_state!`
- `ZernikeBasis`, `compute_zernike!`, `noll_to_nm`
- `Detector`, `DetectorParams`, `DetectorState`, `DetectorExportMetadata`
- `capture!`, `output_frame`, `detector_export_metadata`, `readout_ready`,
  `reset_integration!`
- `DeformableMirror`, `DeformableMirrorParams`, `DeformableMirrorState`,
  `build_influence_functions!`, `apply!`
- `Misregistration`, `apply_misregistration`
- `OPDMap`
- `NCPA`, `NCPABasis`, `KLBasis`, `ZernikeModalBasis`, `M2CBasis`,
  `default_ncpa_basis`
- `SpatialFilter`, `SpatialFilterShape`, `CircularFilter`, `SquareFilter`,
  `FoucaultFilter`, `set_spatial_filter!`, `filter!`

## Atmosphere

- `KolmogorovAtmosphere`, `KolmogorovParams`, `KolmogorovState`
- `MultiLayerAtmosphere`, `MultiLayerParams`, `MultiLayerState`
- `advance!`, `propagate!`
- `update_psd!`, `ensure_psd!`, `phase_screen_von_karman!`
- `phase_variance`, `phase_covariance`, `phase_spectrum`, `covariance_matrix`
- `ft_phase_screen`, `ft_sh_phase_screen`, `PhaseStatsWorkspace`
- `SubharmonicMode`, `FastSubharmonics`, `FidelitySubharmonics`,
  `default_subharmonic_mode`
- `fast_atmosphere`

## Wavefront sensing

- `ShackHartmann`, `ShackHartmannParams`, `ShackHartmannState`
- `PyramidWFS`, `PyramidParams`, `PyramidState`
- `BioEdgeWFS`, `BioEdgeParams`, `BioEdgeState`
- `update_valid_mask!`, `measure!`, `apply_shift_wfs!`, `set_optical_gain!`
- `LiFT`, `lift_interaction_matrix`, `lift_interaction_matrix!`
- `LiFTDampingMode`, `LiFTDampingNone`, `LiFTLevenbergMarquardt`,
  `LiFTAdaptiveLevenbergMarquardt`
- `lgs_elongation_factor`

## Calibration and reconstruction

- `InteractionMatrix`, `interaction_matrix`
- `CalibrationVault`, `with_truncation`
- `ModalBasis`, `KLBasisMethod`, `KLDMModes`, `KLHHtPSD`
- `dm_basis`, `kl_modal_basis`, `modal_basis`, `basis_from_m2c`,
  `basis_projector`
- `AOCalibration`, `ao_calibration`
- `fitting_error`, `fitting_error_dm`
- `GainSensingCamera`, `calibrate!`, `reset_calibration!`,
  `compute_optical_gains!`
- `MetaSensitivity`, `compute_meta_sensitivity_matrix`,
  `estimate_misregistration`, `SPRINT`, `estimate!`
- `ModalReconstructor`, `reconstruct!`, `reconstruct`
- `TomographyNoiseModel`, `RelativeSignalNoise`, `ScalarMeasurementNoise`,
  `DiagonalMeasurementNoise`

## Control

- `AbstractController`
- `DiscreteIntegratorController`, `update!`
- `ClosedLoopRuntime`, `RTCBoundary`, `MultiRTCBoundary`
- `rtc_command`, `rtc_slopes`, `rtc_wfs_frame`, `rtc_science_frame`
- `rtc_wfs_metadata`, `rtc_science_metadata`

## Traits and interfaces

- Abstract interfaces: `AbstractOpticalElement`, `AbstractSource`,
  `AbstractAtmosphere`, `AbstractWFS`, `AbstractDetector`,
  `AbstractDeformableMirror`
- Sensing-mode traits: `SensingMode`, `Diffractive`, `Geometric`
- Detector-noise traits: `NoiseModel`, `NoiseNone`, `NoisePhoton`,
  `NoiseReadout`, `NoisePhotonReadout`
- DM-application traits: `DMApplyMode`, `DMAdditive`, `DMReplace`

## Notes on backend-aware APIs

- Constructors typically accept `T=` and `backend=` keywords.
- Methods accept `AbstractArray` inputs where practical, but structs keep
  concrete array fields so the compiler can specialize.
- FFT-heavy paths use `AbstractFFTs.jl`; future GPU-specialized kernels should
  be introduced behind trait-selected backends rather than by branching on
  symbols in hot loops.

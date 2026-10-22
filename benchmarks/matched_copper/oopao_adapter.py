"""Pinned OOPAO runner for the matched Copper prescription.

Set OOPAO_BACKEND, OOPAO_PRECISION=32, and OOPAO_GPU_RESIDENT before importing
this module or any OOPAO package module. ``upload_opd`` stages float32 OPD in
the selected backend outside a timed step. ``step`` applies that OPD, forms
phase as part of OOPAO propagation, and returns a borrowed native-backend
64×64 photon-rate image. Consume or copy it before the next step.
"""

import os

import numpy as np


def validate_contract(contract):
    """Check the common physical and sampling contract used by both runners."""
    expected = dict(
        resolution=240,
        modulation_points=96,
        modulation_radius_lambda_over_d=5.0,
        diameter_m=1.22,
        wavelength_m=550e-9,
        pupil_samples=30,
        separation_pixels=2,
        edge_pixels=1,
    )
    for key, value in expected.items():
        if contract[key] != value:
            raise ValueError(f"Matched Copper requires {key}={value!r}, got {contract[key]!r}")
    irradiance = contract["photon_irradiance"]
    if not np.isfinite(irradiance) or irradiance <= 0:
        raise ValueError("photon_irradiance must be finite and positive")


class _SingleFFTWorkerPyramid:
    """Mixin that keeps OOPAO's FFT kernel while fixing scipy FFT workers."""

    def _pyramid_transform_block(self, fields, workers=-1, per_field_focal_plane=False):
        # scipy.fft.set_workers does not override this explicit workers=-1.
        # Keep OOPAO's transform implementation and only fix its thread count.
        return super()._pyramid_transform_block(
            fields, workers=1, per_field_focal_plane=per_field_focal_plane
        )


class OopaoRunner:
    """Run OOPAO's Pyramid on a common 240-pixel input and return 64² rates.

    OPD input and resident storage are float32 meters. Pupil arrays use
    ``[row, column]`` ordering and must be binary. The output is summed photons
    per second per binned detector pixel. The OPD copy is staged by
    ``upload_opd``; assigning it to ``Source.OPD_no_pupil`` and forming the
    wavelength-scaled phase happen in ``step``.
    """

    def __init__(self, contract, amplitude, opd, backend="cpu", cpu_workers=1):
        validate_contract(contract)
        if backend not in {"cpu", "cuda"}:
            raise ValueError("backend must be 'cpu' or 'cuda'")
        if not isinstance(cpu_workers, int) or cpu_workers < 1:
            raise ValueError("cpu_workers must be a positive integer")
        if os.environ.get("OOPAO_BACKEND", "").lower() != backend:
            raise RuntimeError(f"Set OOPAO_BACKEND={backend} before importing OOPAO")
        if os.environ.get("OOPAO_PRECISION") != "32":
            raise RuntimeError("Set OOPAO_PRECISION=32 before importing OOPAO")
        expected_resident = "1" if backend == "cuda" else "0"
        if os.environ.get("OOPAO_GPU_RESIDENT") != expected_resident:
            raise RuntimeError(f"Set OOPAO_GPU_RESIDENT={expected_resident} before importing OOPAO")

        amplitude = np.asarray(amplitude, dtype=np.float32)
        if amplitude.shape != (240, 240):
            raise ValueError("Shared amplitude must have shape (240,240)")
        if not np.all((amplitude == 0) | (amplitude == 1)) or not np.any(amplitude):
            raise ValueError("OOPAO custom pupil must be a nonempty binary (240,240) mask")

        # Runtime selection and each OOPAO module's captured backend are
        # process-global, so validate only after the root driver configured it.
        from OOPAO.runtime import array_backend, gpu_resident, precision_bits

        self.xp, gpu_enabled = array_backend()
        if bool(gpu_enabled) != (backend == "cuda"):
            raise RuntimeError("OOPAO's selected array backend does not match the requested backend")
        if gpu_resident() != (backend == "cuda"):
            raise RuntimeError("OOPAO's selected GPU residency does not match the requested backend")
        if precision_bits() != 32:
            raise RuntimeError("OOPAO did not select float32 precision")

        from OOPAO.Source import Source
        from OOPAO.Telescope import Telescope
        from OOPAO.Pyramid import Pyramid

        if backend == "cpu":
            class SingleFFTWorkerPyramid(_SingleFFTWorkerPyramid, Pyramid):
                def initialization(self):
                    # Pyramid assigns its machine-dependent nJobs before calling
                    # initialization; set the runner's requested count before
                    # calibration propagation also runs inside its constructor.
                    self.nJobs = cpu_workers
                    return super().initialization()

            pyramid_type = SingleFFTWorkerPyramid
        else:
            pyramid_type = Pyramid

        self.contract = dict(contract)
        self.backend = backend
        self.cpu_workers = cpu_workers
        self.dtype = np.dtype(np.float32)
        self.dt = np.float32(0.002)
        self.source = Source("V", magnitude=0, display_properties=False)
        # Source's V band is 550 nm in this pinned revision; set explicitly so
        # Pyramid calibration and later phase formation use the contract value.
        self.source.wavelength = np.float32(contract["wavelength_m"])
        self.source.nPhoton = np.float32(contract["photon_irradiance"])
        self.telescope = Telescope(
            resolution=240,
            diameter=1.22,
            samplingTime=float(self.dt),
            pupil=amplitude,
        )
        self.source * self.telescope
        self.sensor = pyramid_type(
            nSubap=30,
            telescope=self.telescope,
            modulation=5.0,
            lightRatio=0,
            postProcessing="fullFrame",
            psfCentering=True,
            n_pix_separation=2,
            calibModulation=50.0,
            n_pix_edge=1,
            nTheta_user_defined=96,
        )
        self.sensor.compute_focal_plane = False
        if backend == "cpu":
            self.sensor.nJobs = cpu_workers
        self.opd_buffer = self.xp.empty((240, 240), dtype=self.xp.float32)
        self.upload_opd(opd)
        self._output_rate = None
        self.tick = 0

        if self.sensor.resolution != 512 or self.sensor.cam.resolution != 64:
            raise RuntimeError(
                "Pinned OOPAO geometry did not produce the required 512 FFT and 64x64 camera"
            )

    def upload_opd(self, opd):
        """Copy a host or backend 240×240 OPD in meters into resident float32 storage."""
        if getattr(opd, "shape", None) != (240, 240):
            raise ValueError("Shared OPD must have shape (240,240)")
        self.opd_buffer[...] = self.xp.asarray(opd, dtype=self.xp.float32)

    def step(self, opd=None):
        """Propagate one resident OPD and return borrowed 64×64 photons/s pixels."""
        if opd is not None:
            self.upload_opd(opd)
        self.source.OPD_no_pupil = self.opd_buffer
        self.sensor.wfs_measure(integrate=False)
        binned_photons = self.sensor.cam.set_binning(
            self.sensor.raw_data, 8, mode="sum"
        )
        self._output_rate = binned_photons / self.dt
        if self._output_rate.shape != (64, 64):
            raise RuntimeError(f"Unexpected OOPAO output shape: {self._output_rate.shape}")
        self.tick += 1
        return self._output_rate

    @property
    def output_rate(self):
        """Borrowed native-backend photon-rate output from the latest step."""
        return self._output_rate

    def metadata(self):
        """Return configuration and allocation-path facts outside timed steps."""
        sensor = self.sensor
        max_batch = sensor._max_batch()
        return dict(
            backend=self.backend,
            precision_bits=32,
            source_wavelength_m=float(self.source.wavelength),
            source_photon_irradiance=float(self.source.nPhoton),
            pupil_shape=tuple(self.telescope.pupil.shape),
            pupil_dtype=str(self.telescope.pupil.dtype),
            pupil_storage="host_numpy",
            source_array_storage="device" if self.source.gpu_resident else "host",
            pyramid_mask_storage=type(sensor.mask).__module__,
            phasor_storage=type(sensor.phasor).__module__,
            fft_dtype=str(sensor.precision_complex),
            fft_side=sensor.resolution,
            output_shape=(64, 64),
            output_units="photons/s/pixel",
            output_binning_factor=8,
            sampling_time_s=float(self.telescope.samplingTime),
            n_jobs=sensor.nJobs,
            cpu_fft_workers=1 if self.backend == "cpu" else None,
            modulation_points=sensor.nTheta,
            modulation_phasors_cached=sensor.modulation_phasors is not None,
            modulation_phasor_shape=(
                tuple(sensor.modulation_phasors.shape)
                if sensor.modulation_phasors is not None
                else None
            ),
            max_batch=max_batch,
            compute_focal_plane=sensor.compute_focal_plane,
            phase_formation="timed_wfs_measure_from_resident_OPD",
            opd_upload="separate_float32_resident_copy",
            stop="OOPAO_native_pyramid_mask",
            centering="four_pixel",
        )

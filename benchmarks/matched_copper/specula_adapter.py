"""Explicit common-prescription boundary for pinned SPECULA ModulatedPyramid.

SPECULA must be initialized before importing this module. This benchmark-only
subclass changes cold optical construction; its propagation and processing-object
lifecycle are inherited unchanged. It is not SPECULA's native prescription.
"""

import numpy as np

from specula.data_objects.electric_field import ElectricField
from specula.data_objects.simul_params import SimulParams
from specula.processing_objects.modulated_pyramid import ModulatedPyramid


class MatchedCopperPyramid(ModulatedPyramid):
    """A common 512 FFT, square stop, outward-prism Copper prescription."""

    def calc_pyr_geometry(self, DpupPix, pixel_pitch, lambda_, FoV,
                          pup_diam, ccd_side, fov_errinf=0.1,
                          fov_errsup=0.5, pup_dist=None, pup_margin=2,
                          fft_res=3.0, min_pup_dist=None):
        # SPECULA's native 10% guard expands the requested FFT. The shared
        # prescription explicitly owns its finite numerical propagation grid.
        if (DpupPix, pup_diam, ccd_side, pup_dist, pup_margin) != (240, 30, 64, 32, 1):
            raise ValueError("Matched Copper requires N=240, diameter=30, distance=32, margin=1, CCD=64")
        if lambda_ != 550.0:
            raise ValueError("Matched Copper requires wavelengthInNm=550.0")
        return dict(fov_res=1, fp_masking=1.0, fft_res=512 / 240,
                    fft_sampling=240, fft_padding=272, fft_totsize=512,
                    tilt_scale=4.0, toccd_side=64, final_ccd_side=64)

    def get_fp_mask(self, totsize, mask_ratio, obsratio=0):
        # Native make_mask always clips a circle, even at ratio=1. The AOS
        # common prescription transmits the complete sampled square plane.
        if totsize != 512 or obsratio != 0:
            raise ValueError("Matched Copper requires a 512 square stop without obstruction")
        return self.xp.ones((512, 512), dtype=self.dtype)

    def get_pyr_tlt(self, p, c):
        # Native folded ramps have +|frequency| phase after removing piston.
        # AOS uses -|frequency|. Negating this cold surface adopts that physical
        # prism orientation, retaining SPECULA's native face construction.
        return -super().get_pyr_tlt(p, c)


def validate_contract(contract):
    expected = dict(resolution=240, modulation_points=96,
                    modulation_radius_lambda_over_d=5.0, diameter_m=1.22,
                    wavelength_m=550e-9, pupil_samples=30,
                    separation_pixels=2, edge_pixels=1)
    for key, value in expected.items():
        if contract[key] != value:
            raise ValueError(f"Matched Copper requires {key}={value!r}, got {contract[key]!r}")
    irradiance = contract["photon_irradiance"]
    if not np.isfinite(irradiance) or irradiance <= 0:
        raise ValueError("photon_irradiance must be finite and positive")


class CopperSpeculaRunner:
    """Own stable input buffers and execute SPECULA's complete optical step.

    Arrays use [row, column] on both sides. OPD is supplied in meters, converted
    to SPECULA phaseInNm in place without a sign change. `step()` returns a
    borrowed resident 64×64 photon-rate array; read/copy it before the next step.
    The native diagnostic can have a different realized internal/output pupil
    sampling and must not be reported as common-prescription equivalence.
    """

    def __init__(self, contract, amplitude, opd, *, target_device_idx=-1,
                 precision=1, matched=True, capture=True):
        validate_contract(contract)
        amplitude = np.asarray(amplitude, dtype=np.float32)
        if amplitude.shape != (240, 240):
            raise ValueError("Shared amplitude must have shape (240,240)")
        if not np.all((amplitude == 0) | (amplitude == 1)) or not np.any(amplitude):
            raise ValueError("Native SPECULA photon normalization requires a nonempty binary pupil")
        self.contract = dict(contract)
        self.matched = matched
        pitch = np.float32(contract["diameter_m"]) / np.float32(contract["resolution"])
        self.ef = ElectricField(240, 240, pitch,
                                S0=np.float32(contract["photon_irradiance"]),
                                target_device_idx=target_device_idx, precision=precision)
        self.xp = self.ef.xp
        self.ef.A[:] = self.xp.asarray(amplitude, dtype=self.ef.dtype)
        self.opd_m = self.xp.empty((240, 240), dtype=self.ef.dtype)
        self.upload_opd(opd)
        simul = SimulParams(pixel_pupil=240, pixel_pitch=pitch)
        sensor_type = MatchedCopperPyramid if matched else ModulatedPyramid
        self.sensor = sensor_type(
            simul_params=simul, wavelengthInNm=550.0,
            fov=(550e-9 / pitch) * (180 * 3600 / np.pi),
            pup_diam=30, output_resolution=64, mod_amp=5.0, mod_step=96,
            mod_type="circular", pup_dist=32, pup_margin=1, fft_res=512 / 240,
            target_device_idx=target_device_idx, precision=precision)
        self.sensor.inputs["in_ef"].set(self.ef)
        self.sensor.stream_enable = capture
        self.sensor.setup()
        self.tick = 0

    def upload_opd(self, opd):
        """Stage physical meters as resident SPECULA nanometers, outside timing."""
        if opd.shape != (240, 240):
            raise ValueError("Shared OPD must have shape (240,240)")
        self.opd_m[:] = self.xp.asarray(opd, dtype=self.ef.dtype)
        self.xp.multiply(self.opd_m, self.ef.dtype(1e9), out=self.ef.phaseInNm)

    def step(self, opd=None):
        """Run input preparation, inherited trigger/capture, and detector rebin."""
        if opd is not None:
            self.upload_opd(opd)
        self.tick += 1
        self.ef.generation_time = self.tick
        if not self.sensor.check_ready(self.tick):
            raise RuntimeError("SPECULA rejected a refreshed resident input")
        self.sensor.trigger()
        self.sensor.post_trigger()
        return self.sensor.out_i.i

    def metadata(self):
        sensor = self.sensor
        probe = np.zeros((2, 2), dtype=np.complex64)
        return dict(
            prescription="common" if self.matched else "native-diagnostic",
            numpy_version=np.__version__,
            numpy_complex64_fft_dtype=str(np.fft.fft2(probe).dtype),
            field_dtype=str(self.ef.field.dtype),
            output_dtype=str(sensor.out_i.i.dtype),
            fft_sampling=sensor.fft_sampling, fft_padding=sensor.fft_padding,
            fft_totsize=sensor.fft_totsize, fft_res=float(sensor.fft_res),
            toccd_side=sensor.toccd_side, final_ccd_side=sensor.final_ccd_side,
            tilt_scale=float(sensor.tilt_scale),
            stop="square" if self.matched else "circular",
            prism_phase="negative_abs_frequency" if self.matched else "positive_abs_frequency",
            cuda_graph= sensor.cuda_graph is not None)

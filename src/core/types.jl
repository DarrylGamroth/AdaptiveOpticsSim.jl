"""Detectors implement capture!(det, psf; rng)."""
abstract type AbstractDetector <: Optics.AbstractOpticalElement end

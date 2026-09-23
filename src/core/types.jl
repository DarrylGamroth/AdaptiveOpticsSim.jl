"""Detectors implement capture!(det, photon_arrival_rate; rng)."""
abstract type AbstractDetector <: Optics.AbstractOpticalElement end

"""Base type for physical optical elements and optical-domain models."""
abstract type AbstractOpticalElement end

"""Telescopes own revisioned pupil geometry used to prepare optical paths."""
abstract type AbstractTelescope <: AbstractOpticalElement end

"""Sources provide wavelength, direction, and declared radiometry."""
abstract type AbstractSource <: AbstractOpticalElement end

"""Application policy for an already formed sampled optical surface."""
abstract type DMApplyMode end

"""Add a sampled surface OPD to the current path OPD."""
struct DMAdditive <: DMApplyMode end

"""Replace the current path OPD with a sampled surface OPD."""
struct DMReplace <: DMApplyMode end

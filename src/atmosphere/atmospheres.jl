"""
    Atmospheres

Canonical owner of atmospheric models and evolving state, source-direction
geometry, prepared directional rendering, and atmosphere-coupled propagation.
"""
module Atmospheres

using KernelAbstractions
using LinearAlgebra
using Random
using SpecialFunctions
using Statistics

import ..AdaptiveOpticsSim:
    AdaptiveOpticsSimError,
    DimensionMismatchError,
    FastProfile,
    FidelityProfile,
    InvalidConfiguration,
    NumericalConditionError,
    ScientificProfile,
    UnsupportedAlgorithm,
    _scaled_kv56_scalar,
    atmosphere_profile,
    default_build_backend,
    default_fidelity_profile,
    fftfreq!,
    materialize_build

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    ExecutionStyle,
    KernelLaunchPhase,
    ScalarCPUStyle,
    _resolve_array_backend,
    _resolve_backend_selector,
    backend,
    begin_kernel_phase,
    compute_device,
    execute_fft_plan!,
    execution_style,
    finish_kernel_phase!,
    launch_kernel!,
    launch_kernel_async!,
    plan_ifft_backend!,
    queue_kernel!,
    randn_backend!,
    require_same_backend,
    synchronize_backend!

import ..Optics:
    AbstractOpticalElement,
    AbstractSource,
    AchromaticSpectralCoordinate,
    Asterism,
    ElectricField,
    ExtendedSource,
    FresnelPropagation,
    MetricCoordinates,
    NonCombinableProduct,
    OpticalPlaneMetadata,
    PointSampledMeasure,
    PupilFieldFormationPlan,
    PupilFunction,
    PupilPlane,
    SpectralSource,
    Telescope,
    UnspecifiedNormalization,
    _accumulate_field_intensity!,
    _accumulate_field_intensity_async!,
    _converted_nonnegative_finite,
    _converted_positive_finite,
    apply_phase_async!,
    coordinates_xy_arcsec,
    extended_source_asterism,
    fill_electric_field_async!,
    freeze_source,
    intensity!,
    plane_metadata,
    prepare_pupil_field,
    propagate_field!,
    pupil_mask,
    pupil_reflectivity,
    require_leaf_source,
    require_same_plane_grid,
    source_height_m,
    source_radiometric_value,
    source_with_wavelength_and_radiometric_value,
    spectral_bundle,
    validate_plane_storage,
    wavelength

"""Invalid or non-monotonic explicit atmosphere model time."""
struct AtmosphereTimeError <: AdaptiveOpticsSimError
    msg::String
end

"""Missing, stale, or incompatible atmosphere epoch identity."""
struct AtmosphereEpochError <: AdaptiveOpticsSimError
    msg::String
end

"""Base type for atmospheric models."""
abstract type AbstractAtmosphere <: AbstractOpticalElement end

"""
Atmospheres whose shared evolution is advanced with explicit model time and
published as immutable current-state epoch identity tokens.
"""
abstract type AbstractTimedAtmosphere <: AbstractAtmosphere end

# Extension seam retained for untimed/static test atmospheres and user models.
# Maintained stochastic atmospheres use `advance_by!` / `advance_to!` instead.
function advance! end

include("source_geometry.jl")
include("kolmogorov.jl")
include("infinite_screen_math.jl")
include("infinite_screen.jl")
include("multilayer.jl")
include("direction_batch.jl")
include("phase_stats.jl")
include("atmospheric_field_propagation.jl")
include("api.jl")

end # module Atmospheres

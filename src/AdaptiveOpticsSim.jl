module AdaptiveOpticsSim

__precompile__(true)

using AbstractFFTs
import FFTW
using KernelAbstractions
using LinearAlgebra
using Logging
using Random
using SparseArrays
using SpecialFunctions
using Statistics

import Base: filter!

"""
AdaptiveOpticsSim.jl

Julia adaptive optics simulation toolkit (in development).
"""
const PROJECT_STATUS = :in_development

include("core/errors.jl")
include("core/profiles.jl")
include("core/random_services.jl")
include("backends/backends.jl")

# The remaining flat source files are migrated owner by owner. Until their
# modules land, import the canonical Backends identities strictly for
# package-internal composition; the root API does not export or mark them
# public.
import .Backends:
    AMDGPUBackend,
    AMDGPUBackendTag,
    AbstractArrayBackend,
    AbstractComputeDevice,
    AbstractReductionExecutionPlan,
    AcceleratorComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    CUDABackend,
    CUDABackendTag,
    DirectReductionPlan,
    ExecutionStyle,
    GPUBackendTag,
    GPUPrecisionPolicy,
    HostComputeDevice,
    HostMirrorReductionPlan,
    KernelLaunchPhase,
    MetalBackend,
    MetalBackendTag,
    ScalarCPUStyle,
    SplitGPUPrecision,
    UnifiedGPUPrecision,
    _AbstractPreparedDeviceExecutionContext,
    _HostPreparedDeviceExecutionContext,
    _clamp_array!,
    _poisson_noise!,
    _prepare_device_execution_context,
    _prepared_device_execution_compute_device,
    _rand_uniform_backend!,
    _randn_backend!,
    _resolve_array_backend,
    _resolve_backend_selector,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    _write_integer_output!,
    allocate_array,
    array_backend_selector,
    array_backend_type,
    available_gpu_backends,
    backend,
    backend_fill,
    backend_matmul,
    backend_matmul_transpose_right,
    backend_maximum_value,
    backend_rand,
    backend_randn,
    backend_sum_value,
    backend_symmetric_product,
    backend_type,
    backend_zeros,
    begin_kernel_phase,
    clamp_array!,
    compute_device,
    compute_device_backend,
    compute_device_identifier,
    default_gpu_precision_policy,
    disable_scalar_backend!,
    execute_fft_plan!,
    execution_style,
    finish_kernel_phase!,
    gpu_backend_array_type,
    gpu_backend_loaded,
    gpu_backend_name,
    gpu_build_type,
    gpu_runtime_type,
    high_accuracy_gpu_precision_policy,
    host_array,
    launch_kernel!,
    launch_kernel_async!,
    masked_sum2d,
    packed_valid_pair_mean,
    plan_fft_backend!,
    plan_ifft_backend!,
    poisson_noise!,
    poisson_noise_async!,
    queue_kernel!,
    rand_uniform_backend!,
    randn_backend!,
    randn_backend_async!,
    reduction_execution_plan,
    reduction_host_view,
    reduction_parent_source,
    require_same_backend,
    resolve_array_backend,
    same_backend,
    set_fft_provider_threads!,
    synchronize_backend!,
    use_host_build_algebra,
    write_integer_output!

include("core/inverse_policies.jl")
include("core/utils.jl")

# Configuration serialization remains a root-owned cross-domain seam. Declare
# its generic before domain modules add methods, then install the shared
# fallbacks after the first physical owner has loaded.
function config_value end

include("optics/optics.jl")

# The remaining flat source files are migrated owner by owner. Routine Optics
# vocabulary is imported for package-internal composition only; these imports
# are not root exports or qualified-public compatibility bindings.
using .Optics
import .Optics:
    AbstractOpticalElement,
    AbstractPropagationModel,
    AbstractSource,
    AbstractSourceCompositionStyle,
    AbstractSpectralCoordinate,
    AbstractTelescope,
    DMApplyMode,
    ExpandedSourceComposition,
    ExtendedSource,
    LeafSourceComposition,
    NCPA,
    PreparedBundledDirectImaging,
    PreparedDirectImaging,
    PreparedIncoherentDirectImaging,
    LGSProfileNaProfile,
    LGSProfileNone,
    _FixedOpticalProductVector,
    _accumulate_field_intensity!,
    _accumulate_field_intensity_async!,
    _converted_nonnegative_finite,
    _converted_positive_finite,
    _direct_imaging_batch_product_contract,
    _pupil_cell_area,
    _pupil_diameter_m,
    _pupil_resolution,
    _require_physical_photon_irradiance,
    aperture_revision,
    apply_centering_phase!,
    apply_phase_async!,
    axis_centering,
    centered_grid_origin,
    coordinates_xy_arcsec,
    electric_field_wavelength,
    field_embedding_offsets,
    freeze_source,
    fraunhofer_intensity_stack!,
    intensity!,
    is_leaf_source,
    is_lgs_source,
    lgs_elongation_factor,
    lgs_profile,
    photon_irradiance,
    pixel_mask_grid,
    plane_metadata,
    require_leaf_source,
    require_metric_coordinates,
    require_same_plane_grid,
    require_centered_plane_geometry,
    source_measurement_signature,
    source_composition_style,
    source_height_m,
    source_with_wavelength_and_radiometric_value,
    spectral_bundle,
    surface_opd,
    wavelength,
    validate_plane_storage

include("core/types.jl")
include("core/config.jl")
include("core/kv56.jl")
include("core/workspace.jl")
include("core/parallel.jl")
include("core/telemetry.jl")

include("optics/zernike.jl")
include("optics/misregistration.jl")
include("optics/controllable_optics.jl")
include("optics/deformable_mirror.jl")
include("detectors/detector.jl")
include("optics/spatial_filter.jl")
include("atmosphere/source_geometry.jl")
include("atmosphere/kolmogorov.jl")
include("atmosphere/infinite_screen_math.jl")
include("atmosphere/infinite_screen.jl")
include("atmosphere/multilayer.jl")
include("atmosphere/direction_batch.jl")
include("atmosphere/phase_stats.jl")
include("optics/atmospheric_field_propagation.jl")
include("calibration/modal_basis.jl")
include("calibration/ncpa.jl")
include("wfs/sensing_modes.jl")
include("wfs/stage_contracts.jl")
include("wfs/grouped.jl")
include("wfs/calibration.jl")
include("wfs/elongation.jl")
include("wfs/subapertures.jl")
include("wfs/focal_plane_modulation.jl")
include("wfs/shack_hartmann.jl")
include("wfs/pyramid.jl")
include("wfs/bioedge.jl")
include("wfs/zernike.jl")
include("wfs/curvature.jl")
include("wfs/interface.jl")
include("wfs/lift.jl")
include("plant/plant.jl")
include("calibration/interaction_matrix.jl")
include("control/controller.jl")
include("control/reconstructors.jl")
include("calibration/control_matrix.jl")
include("calibration/fitting_error.jl")
include("calibration/ao_calibration.jl")
include("calibration/gain_sensing_camera.jl")
include("calibration/misregistration_identification.jl")
include("calibration/ad_sensitivities.jl")
include("control/ensemble.jl")
include("tomography/parameters.jl")
include("tomography/geometry.jl")
include("tomography/fitting.jl")
include("tomography/reconstructors.jl")

include("exports.jl")

end # module AdaptiveOpticsSim

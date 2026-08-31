"""
    Optics

Canonical owner of foundational optical locations, coordinates, products,
sources, telescope/aperture models, backend-portable propagation, explicit
sampled NCPA, controllable optics, and reusable physical WFS components.
"""
module Optics

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault
using KernelAbstractions
using LinearAlgebra
using Statistics

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration,
    UnsupportedAlgorithm,
    _fixed_size_union_vector,
    _union_members_are_concrete,
    circshift2d!,
    config_value,
    fftfreq!,
    fftshift2d!

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    ExecutionStyle,
    ScalarCPUStyle,
    _throw_compute_device_error,
    _with_compute_device,
    _resolve_array_backend,
    _resolve_backend_selector,
    allocate_array,
    backend,
    begin_kernel_phase,
    compute_device,
    compute_device_backend,
    execute_fft_plan!,
    execution_style,
    finish_kernel_phase!,
    launch_kernel!,
    launch_kernel_async!,
    plan_fft_backend!,
    plan_ifft_backend!,
    plan_repeated_fft_backend!,
    plan_repeated_ifft_backend!,
    queue_kernel!,
    require_same_backend,
    synchronize_backend!

include("types.jl")
include("aperture_masks.jl")
include("telescope.jl")
include("source.jl")
include("spectrum.jl")
include("planes.jl")
include("electric_field.jl")
include("propagation.jl")
include("asterism.jl")
include("extended_source.jl")
include("direct_imaging.jl")
include("direct_imaging_batch.jl")
include("zernike.jl")
include("opd_map.jl")
include("misregistration.jl")
include("ncpa.jl")
include("controllable_optics.jl")
include("deformable_mirror.jl")

# Canonical Optics-owned operation; do not extend Base.filter! because this
# contract transforms optical products rather than filtering a collection.
function filter! end

include("spatial_filter.jl")
include("focal_plane_modulation.jl")
include("microlens_array.jl")
include("wfs_components.jl")
include("api.jl")

end # module Optics

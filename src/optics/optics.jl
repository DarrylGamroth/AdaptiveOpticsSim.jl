"""
    Optics

Canonical owner of foundational optical locations, coordinates, products,
sources, telescope/aperture models, backend-portable propagation, and explicit
sampled NCPA. Reusable controllable and WFS optical components complete their
move in the following namespace gate.
"""
module Optics

using KernelAbstractions
using Statistics

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration,
    UnsupportedAlgorithm,
    circshift2d!,
    config_value,
    fftfreq!

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorStyle,
    CPUBackend,
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
    plan_fft_backend!,
    plan_ifft_backend!,
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
include("opd_map.jl")
include("ncpa.jl")
include("api.jl")

end # module Optics

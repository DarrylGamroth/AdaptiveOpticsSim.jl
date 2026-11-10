"""
    Optics

Canonical owner of foundational optical locations, coordinates, products,
sources, telescope/aperture models, backend-portable propagation, explicit
sampled NCPA, controllable optics, and reusable physical WFS components.
"""
module Optics

# Optional Proper implementation; data readers are separate optional extensions.
"""
    prepare_compiled_optical_paths(metadata, read_array; target, T=Float64,
        dm_opd=nothing, entrance_opd=nothing, entrance_attenuation=nothing,
        summation=:path_pairs)

Prepare a frozen multi-path Fresnel and spectral model. Load Proper.jl to
activate execution. `read_array(reference)` supplies canonical logical[y,x]
matrices during cold preparation; execution retains no metadata dictionary or
reader. OPD inputs are metres; entrance attenuation is the finite additional
log-amplitude quantity `a` in `exp(-a)` on the compiled native input grid.
Backend execution policies, including CUDA math mode, are fixed at preparation
and capture; changing them requires a new preparation/capture.
"""
function prepare_compiled_optical_paths end
"""
    propagate_compiled_optical_paths!(owner)

Update a standalone compiled optical model's fringed, unfringed and reference
relative-intensity products using its bound native-grid live inputs. Returns the
retained products after completion. Exposure timing, detector acquisition and
physical radiometry remain the caller's responsibility. Requires Proper.jl.
"""
function propagate_compiled_optical_paths! end

"""
    fitted_scc_filter_graph(model_path, model_sha256)

Construct an opt-in FilterGraphs definition and its Algorithm declarations for
a fitted SCC optical model and a downstream pixel MTF. Inputs are the native
DM and entrance OPD in metres and the entrance log-amplitude attenuation.
Load Proper, FITSIO, YAML and FilterGraphAlgorithms before calling. Prepare
the returned dictionary with FilterGraphs and select native capture on GPU.
"""
function fitted_scc_filter_graph end

"""
    llowfs_filter_graph(definition)

Build a FilterGraphs description and its Algorithm declarations from a Float32
AOS LLOWFS optical train followed by pixel integration. The graph accepts
real pupil OPD, pupil amplitude and DM OPD, and exports real detector-grid
intensity. Load Proper, FilterGraphs and FilterGraphAlgorithms before calling.
"""
function llowfs_filter_graph end

"""
    pwfs_filter_graph(definition)

Build a FilterGraphs description and Algorithm declarations for a Float32
Pyramid photon-rate node with real pupil OPD input. The prepared Algorithm owns
its diffraction and modulation workspaces. Load Proper, FilterGraphs and
FilterGraphAlgorithms before calling.
"""
function pwfs_filter_graph end

"""
    with_filtergraphs_pwfs(definition; node=:pyramid)

Replace the named Float32 OPD-input Pyramid node with the FilterGraphAlgorithms
declaration, retaining the graph's links and boundaries. Prepare and step the
result with AOS to retain source, exposure, detector, and DM timing. On GPU the
node is prepared against the AOS graph's device buffers and records into the
enclosing native capture.
"""
function with_filtergraphs_pwfs end

"""
    prepare_filtergraphs_external_pupil(prepared; single_scc_product=false,
        merge_scc_fields=false, llowfs_definition=nothing)

Select the experimental FilterGraphs optical executor for an unused, fitted
SPIDERS SCC/LLOWFS external pupil. The returned owner takes over the original
optical buffers, exposure history and detector response. Retain and step only
the returned owner. Keep using the existing external pupil and DM input API.

Load FilterGraphs, FilterGraphAlgorithms and Proper to enable this method.
CPU uses a prepared FilterGraphs graph. Supported accelerators require native
capture of each SCC chopping selection. The caller and the external pupil own
the exposure and feedback clocks.

Set `single_scc_product=true` to select the experimental CUDA path for a
fitted SCC observation with one compiled optical node and no pixel MTF.
It computes only the intensity used by the current exposure. Combining the
complex fields before their common transforms changes Float32 roundoff.
Set `merge_scc_fields=true` with `single_scc_product=true` to merge the fringed
fields after the second differing mask and use only the closed field for the
unfringed image when the compiled prescription admits that traversal.

Supply the authored LLOWFS graph definition as `llowfs_definition` to prepare
its two direct Algorithms once in their own FilterGraphs graph. Use the same
prescription as the prepared AOS observation. The returned owner feeds the
detector-grid intensity to the existing exposure accumulator.
"""
function prepare_filtergraphs_external_pupil end

# Backend policy hooks for cold preparation; other backends retain the ordinary
# multiplier storage and BLAS arithmetic.
_prepared_optical_train_enabled(target, field) = false
_compiled_paths_branch_batch_enabled(target, field) = false
_compiled_paths_resample_batch_size(target, field, image_size, count) = 1
function _compiled_paths_resample_mul!(field, temporary, kernel, coefficients)
    matrix = view(kernel, :, :, 1)
    for i in axes(field, 3)
        input = view(field, :, :, i)
        work = view(temporary, :, :, i)
        _enqueue_mul!(work, matrix, input, coefficients)
        _enqueue_mul!(input, work, transpose(matrix), coefficients)
    end
    return nothing
end
_compiled_paths_edge_cache_enabled(target, field) = false
_prepare_compiled_paths_mul_coefficients(field, coefficients) = coefficients

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault, FixedSizeMatrixDefault
using KernelAbstractions
using LinearAlgebra
using Statistics
import Interpolations

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
    _enqueue_mul!,
    _throw_compute_device_error,
    _with_compute_device,
    _resolve_array_backend,
    _resolve_backend_selector,
    allocate_array,
    backend,
    begin_kernel_phase,
    compute_device,
    compute_device_backend,
    copyto_backend_async!,
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

include("zemax.jl")
include("types.jl")
include("aperture_masks.jl")
include("telescope.jl")
include("source.jl")
include("spectrum.jl")
include("spectral_synthesis.jl")
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
include("second_order_response.jl")
include("deformable_mirror.jl")
include("measured_influences.jl")

# Canonical Optics-owned operation; do not extend Base.filter! because this
# contract transforms optical products rather than filtering a collection.
function filter! end

include("spatial_filter.jl")
include("focal_plane_modulation.jl")
include("microlens_array.jl")
include("wfs_components.jl")
include("api.jl")

end # module Optics

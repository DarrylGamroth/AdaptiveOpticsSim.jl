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
    shack_hartmann_rate_filter_graph(definition)

Return a FilterGraphs recipe for an AOS Float32 Shack–Hartmann optical-rate
node. The graph accepts pupil OPD and produces the photon-rate mosaic.
"""
function shack_hartmann_rate_filter_graph end

"""
    cmos_filter_graph(definition)

Return a FilterGraphs recipe for an AOS Float32 CMOS acquisition node.
The graph accepts photon rate and produces the detector frame using the
configured exposure and seeded noise model.
"""
function cmos_filter_graph end

"""
    shack_hartmann_camera_filter_graph(definition; frame_output=:frame)

Build the Float32 Shack–Hartmann photon-rate, CMOS acquisition and UInt16
frame-layout camera that produces the selected AOS output. The standalone
FilterGraphs recipe accepts pupil OPD and publishes a complete raw frame.
The CMOS uses its configured exposure. Supply the corresponding timestamp
and exposure in input sample metadata.
"""
function shack_hartmann_camera_filter_graph end

"""
    with_filtergraphs_shack_hartmann_camera(definition; frame_output=:frame)

Replace the selected Shack–Hartmann camera with FilterGraphAlgorithms
nodes at their original AOS graph positions. Retain surrounding model-time
execution, input and diagnostic routes, and DM command timing. Prepare
with native capture on GPU.
"""
function with_filtergraphs_shack_hartmann_camera end

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
    pwfs_camera_filter_graph(definition)

Build a FilterGraphs description for the real Float32 OPD-input Pyramid,
EMCCD and UInt16 frame layout that produce the AOS `:frame` output.
The EMCCD integrates for the exposure duration in its AOS configuration.
The graph propagates sample metadata supplied by the caller.
"""
function pwfs_camera_filter_graph end

"""
    with_filtergraphs_pwfs_camera(definition)

Bind the FilterGraphAlgorithms camera nodes at their original positions in
the AOS graph. Prepare with the existing AOS API and select native capture
on GPU. AOS retains the source, model clock, exposure and DM command
boundaries.
"""
function with_filtergraphs_pwfs_camera end

"""
    gpi2_plant_filter_graph(definition)

Build a FilterGraphs recipe from a fitted GPI2 plant. The recipe retains
the authored direct routes and diagnostic outputs. The science relay,
atmosphere, mirror surfaces, pupil composition, and PWFS camera use the
AOS physical models. Sparse calibration arrays are captured at preparation.
"""
function gpi2_plant_filter_graph end

"""
    prepare_filtergraphs_gpi2_hil(definition; command_input, target, exposure_ns)

Prepare the GPI2 plant on FilterGraphs with the AOS lockstep HIL interface.
The caller owns the model-time driver and stages commands through
`hil_command_buffer`. A completed sample requires the matching command
response before the next frame. The optional `exposure_ns` must match the
authored detector exposure. Accelerator targets require native capture.
"""
function prepare_filtergraphs_gpi2_hil end

"""
    prepare_filtergraphs_pwfs_hil(definition; command_input, frame_slots=nothing)

Prepare the authored CPU simulator with its PWFS camera executed by
FilterGraphs. Supply a tuple of caller-owned raw frame arrays. Step with
`step_hil_frame_at!(owner, driver, slot)`, selecting an idle slot, and use
`hil_frame_buffer(owner)` after a successful step. The command buffer and
sequence adoption follow the AOS lockstep HIL contract. AOS retains the
source and model clock. Use a FilterGraphs build with CPU loan pools.

Omit `frame_slots` when the caller obtains a new output loan each step. Pass
the idle loan to `step_hil_frame_at!` as a `Matrix{UInt16}` with the reversed
logical shape. The camera writes the logical frame in row-major wire order.
The destination must be independent of all simulator storage and remain valid
for the call. In this mode `hil_frame_buffer(owner)` describes the frame
format only.
"""
function prepare_filtergraphs_pwfs_hil end

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
using Base.ScopedValues: ScopedValue, with

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

const _propagation_checkpoint = ScopedValue{Union{Nothing,Function}}(nothing)

"""
    with_propagation_checkpoint(f, checkpoint)

Run `f()` with a cooperative callback between CPU compiled optical propagation
operations on the current task. The callback may pause or throw before the
step publishes its products. After a throw, reset a prepared graph before
reusing it. The callback must not access or mutate the in-flight optics.
Accelerator execution does not invoke this callback.
"""
function with_propagation_checkpoint(f::F, checkpoint::C) where {F,C<:Function}
    return with(f, _propagation_checkpoint => checkpoint)
end

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

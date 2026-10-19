function _pixel_mtf_shape(shape::Tuple{Integer,Integer})
    all(n -> !(n isa Bool) && 0 < n <= typemax(Int), shape) ||
        throw(InvalidConfiguration("pixel MTF shape must contain two positive integer dimensions"))
    return (Int(shape[1]), Int(shape[2]))
end

_pixel_mtf_shape(shape) = throw(InvalidConfiguration(
    "pixel MTF shape must be a tuple of two positive integer dimensions"))

function _pixel_mtf_nonnegative_parameter(value::Real, name)
    isfinite(value) && value >= 0 || throw(InvalidConfiguration(
        "$name must be finite and nonnegative"))
    converted = Float64(value)
    isfinite(converted) || throw(InvalidConfiguration(
        "$name must remain finite at Float64 precision"))
    return converted
end

_pixel_mtf_nonnegative_parameter(value, name) = throw(InvalidConfiguration(
    "$name must be a real finite nonnegative value"))

"""
    pixel_mtf_transfer(shape; box_width_pixels, diffusion_variance_pixels2,
                       ipc_per_neighbor, T=Float64)

Build a host `Matrix{T}` on the declared two-dimensional frequency grid in
corner-origin FFT ordering, with frequencies in cycles per pixel and unit DC.
`T` must be `Float32` or `Float64`. All three response parameters are explicit:
the square integration width in pixel pitches, Gaussian diffusion variance in
pixel pitches squared, and the fraction coupled to each of four neighbours.
Width and variance must be finite and nonnegative; coupling must lie in [0, 0.25].

For row/column frequencies `(fy, fx)`, the real transfer is
`sinc(w*fy)*sinc(w*fx)*exp(-2π²*v*(fy²+fx²)) *
(1-4a+2a*(cos(2π*fy)+cos(2π*fx)))`. Fit parameters describe an effective curve:
its component contributions may be degenerate and do not independently measure
integration, diffusion, or charge coupling.

Pass the result to [`prepare_pixel_mtf`](@ref) with the same padded shape and
pixel-coordinate convention. Choose the box term to account for any integration
already present in optical sampling or detector acquisition. This builder sets
no detector defaults and supplies no noise covariance. It does not guarantee a
nonnegative finite-grid inverse kernel; application does not clamp the output.
"""
function pixel_mtf_transfer(shape; box_width_pixels, diffusion_variance_pixels2,
        ipc_per_neighbor, T=Float64)
    T in (Float32, Float64) || throw(InvalidConfiguration(
        "pixel MTF transfer type must be Float32 or Float64"))
    nrows, ncols = _pixel_mtf_shape(shape)
    width = _pixel_mtf_nonnegative_parameter(box_width_pixels, "box_width_pixels")
    variance = _pixel_mtf_nonnegative_parameter(diffusion_variance_pixels2,
        "diffusion_variance_pixels2")
    coupling = _pixel_mtf_nonnegative_parameter(ipc_per_neighbor, "ipc_per_neighbor")
    coupling <= 0.25 || throw(InvalidConfiguration(
        "ipc_per_neighbor must lie in [0, 0.25]"))
    # Even grids place the negative Nyquist frequency at index n/2 + 1;
    # odd grids retain their last positive frequency before the negative branch.
    fy = [(k <= (nrows - 1) ÷ 2 ? k : k - nrows) / nrows for k in 0:(nrows - 1)]
    fx = [(k <= (ncols - 1) ÷ 2 ? k : k - ncols) / ncols for k in 0:(ncols - 1)]
    transfer = Matrix{T}(undef, nrows, ncols)
    for column in 1:ncols, row in 1:nrows
        y, x = fy[row], fx[column]
        box = sinc(width * y) * sinc(width * x)
        diffusion = exp(-2pi^2 * (variance * (y^2 + x^2)))
        capacitance = 1 - 4coupling + 2coupling * (cospi(2y) + cospi(2x))
        transfer[row, column] = T(box * diffusion * capacitance)
    end
    return transfer
end

"""
    PixelMTFPlan

Run-immutable grid, padding, target, and owned calibration for a pixel MTF
filter. The transfer uses corner-origin FFT frequency ordering and unit DC.
FFT execution resources and scratch belong to the workspace.
"""
struct PixelMTFPlan{A,B,D,S}
    input_shape::Tuple{Int,Int}
    padded_shape::Tuple{Int,Int}
    padding::Int
    transfer::A
    backend::B
    device::D
    style::S
end

struct PixelMTFWorkspace{A,F,I}
    buffer::A
    fft_plan::F
    ifft_plan::I
end

struct PixelMTFProducts{A}
    output::A
end

struct _PixelMTFPreparationToken end
const _PIXEL_MTF_PREPARATION_TOKEN = _PixelMTFPreparationToken()

"""
    PreparedPixelMTF

Single-writer binding of exact input and output storage, an owned pixel MTF
calibration, and replaceable FFT workspace. Construct with `prepare_pixel_mtf`.
Input samples may change between calls; shape and storage bindings remain fixed.
"""
struct PreparedPixelMTF{I,P,W,O}
    input::I
    plan::P
    workspace::W
    products::O

    function PreparedPixelMTF(::_PixelMTFPreparationToken,
            input::I, plan::P, workspace::W, products::O) where {I,P,W,O}
        return new{I,P,W,O}(input, plan, workspace, products)
    end
end

"""
    prepare_pixel_mtf(input::AbstractMatrix, transfer::AbstractMatrix;
                      padding=16, output=similar(input))

Prepare a deterministic pixel MTF filter on `Float32` or `Float64` intensity
storage. `transfer` is a real finite matrix of size `size(input) .+ 2padding`,
in corner-origin FFT frequency order, with its DC value `transfer[1, 1]` equal
to one within eight machine epsilons of the input precision. It is copied to
the input's compute device during preparation; host calibration is accepted.

The equally typed input and output must have the same shape and exact compute
device, one-based indexing, and no aliases with each other or the calibration.
The operation zero-pads, applies the transfer by FFT, and crops the real inverse
transform. It does not impose a noise model or modify detector defaults.

Preparation may allocate and synchronize. Repeated `apply_pixel_mtf!` calls
reuse workspace and enqueue accelerator work on the caller's current stream
on the prepared device;
the caller owns completion before reading output or changing bound storage.
"""
function prepare_pixel_mtf(input::AbstractMatrix, transfer::AbstractMatrix;
        padding=16, output=similar(input))
    T = eltype(input)
    T in (Float32, Float64) || throw(InvalidConfiguration(
        "pixel MTF input must have Float32 or Float64 elements"))
    output isa AbstractMatrix && eltype(output) === T || throw(InvalidConfiguration(
        "pixel MTF output must be a matrix with the input element type"))
    eltype(transfer) <: Real || throw(InvalidConfiguration(
        "pixel MTF transfer must have real elements"))
    padding isa Integer && !(padding isa Bool) && 0 <= padding <= typemax(Int) ÷ 2 ||
        throw(InvalidConfiguration("pixel MTF padding must be a nonnegative integer"))
    pad = Int(padding)
    Base.require_one_based_indexing(input, output, transfer)
    shape = size(input)
    all(>(0), shape) || throw(InvalidConfiguration("pixel MTF input must be nonempty"))
    all(n -> n <= typemax(Int) - 2pad, shape) || throw(InvalidConfiguration(
        "pixel MTF padded dimensions overflow Int"))
    padded_shape = (shape[1] + 2pad, shape[2] + 2pad)
    size(output) == shape || throw(DimensionMismatchError(
        "pixel MTF output shape must match input"))
    size(transfer) == padded_shape || throw(DimensionMismatchError(
        "pixel MTF transfer shape must match the padded input"))
    (Base.mightalias(input, output) || Base.mightalias(input, transfer) ||
        Base.mightalias(output, transfer)) && throw(InvalidConfiguration(
        "pixel MTF input, output, and supplied transfer must not alias"))
    require_same_backend(input, output)
    device = compute_device(input)
    output_device = compute_device(output)
    output_device == device || _throw_compute_device_error(
        :prepare_pixel_mtf, :wrong_device, device,
        "pixel MTF output occupies $(output_device)")
    # Validate both the supplied values and their representation on the target.
    # Conversion can overflow even when the original Float64 values are finite.
    host_transfer = T.(host_array(transfer))
    all(isfinite, host_transfer) || throw(InvalidConfiguration(
        "pixel MTF transfer must remain finite at input precision"))
    isapprox(host_transfer[1, 1], one(T); rtol=8eps(T), atol=zero(T)) ||
        throw(InvalidConfiguration("pixel MTF transfer must have unit DC"))
    owned_transfer, workspace = _with_compute_device(device) do
        owned_transfer = allocate_device_array(device, T, padded_shape...)
        copyto!(owned_transfer, host_transfer)
        buffer = allocate_device_array(device, Complex{T}, padded_shape...)
        fft_plan = plan_fft_backend!(buffer)
        ifft_plan = plan_ifft_backend!(buffer)
        return owned_transfer, PixelMTFWorkspace(buffer, fft_plan, ifft_plan)
    end
    plan = PixelMTFPlan(shape, padded_shape, pad, owned_transfer,
        backend(input), device, execution_style(input))
    return PreparedPixelMTF(_PIXEL_MTF_PREPARATION_TOKEN, input, plan,
        workspace, PixelMTFProducts(output))
end

@kernel function _pixel_mtf_pad_kernel!(buffer, input, padding, nrows, ncols)
    i, j = @index(Global, NTuple)
    if padding < i <= padding + nrows && padding < j <= padding + ncols
        @inbounds buffer[i, j] = input[i - padding, j - padding]
    else
        @inbounds buffer[i, j] = zero(eltype(buffer))
    end
end

@kernel function _pixel_mtf_multiply_kernel!(buffer, transfer)
    i, j = @index(Global, NTuple)
    @inbounds buffer[i, j] *= transfer[i, j]
end

@kernel function _pixel_mtf_crop_kernel!(output, buffer, padding)
    i, j = @index(Global, NTuple)
    @inbounds output[i, j] = real(buffer[i + padding, j + padding])
end

function _pixel_mtf_pad!(::ScalarCPUStyle, buffer, input, plan)
    fill!(buffer, zero(eltype(buffer)))
    pad = plan.padding
    # Preparation validates one-based axes and every accessed shape.
    @inbounds for j in axes(input, 2), i in axes(input, 1)
        buffer[i + pad, j + pad] = input[i, j]
    end
    return buffer
end

function _pixel_mtf_pad!(style::AcceleratorStyle, buffer, input, plan)
    launch_kernel_async!(style, _pixel_mtf_pad_kernel!, buffer, input,
        plan.padding, plan.input_shape...; ndrange=plan.padded_shape)
    return buffer
end

function _pixel_mtf_multiply!(::ScalarCPUStyle, buffer, transfer, plan)
    @inbounds for j in axes(buffer, 2), i in axes(buffer, 1)
        buffer[i, j] *= transfer[i, j]
    end
    return buffer
end

function _pixel_mtf_multiply!(style::AcceleratorStyle, buffer, transfer, plan)
    launch_kernel_async!(style, _pixel_mtf_multiply_kernel!, buffer, transfer;
        ndrange=plan.padded_shape)
    return buffer
end

function _pixel_mtf_crop!(::ScalarCPUStyle, output, buffer, plan)
    pad = plan.padding
    @inbounds for j in axes(output, 2), i in axes(output, 1)
        output[i, j] = real(buffer[i + pad, j + pad])
    end
    return output
end

function _pixel_mtf_crop!(style::AcceleratorStyle, output, buffer, plan)
    launch_kernel_async!(style, _pixel_mtf_crop_kernel!, output, buffer,
        plan.padding; ndrange=plan.input_shape)
    return output
end

"""
    apply_pixel_mtf!(owner::PreparedPixelMTF)

Filter the bound intensity into `owner.products.output`, returning that exact
output matrix. Accelerator execution is asynchronous on the caller's current
stream on the prepared device; the caller must establish completion before
reading or changing data.
The warmed CPU execution reuses all storage and FFT plans.
"""
function apply_pixel_mtf!(owner::PreparedPixelMTF)
    plan = owner.plan
    workspace = owner.workspace
    buffer = workspace.buffer
    _pixel_mtf_pad!(plan.style, buffer, owner.input, plan)
    enqueue_fft_plan!(buffer, workspace.fft_plan)
    _pixel_mtf_multiply!(plan.style, buffer, plan.transfer, plan)
    enqueue_fft_plan!(buffer, workspace.ifft_plan)
    return _pixel_mtf_crop!(plan.style, owner.products.output, buffer, plan)
end

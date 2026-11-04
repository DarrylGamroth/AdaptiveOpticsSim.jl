# Prepared normalization launches use 256 lanes, cld(length(field), 256)
# partials and one power element. The field guard and fixed reduction tree
# bound every access, including padded lanes and local_i + stride.
# Two fixed reduction stages retain every result on the device. Each workgroup
# reduces in one defined tree order; no atomics or device-to-host scalar queries.
KernelAbstractions.@kernel function _compiled_paths_power_partials_kernel!(partials, field, count)
    i = @index(Global, Linear)
    local_i = @index(Local, Linear)
    group = @index(Group, Linear)
    shared = @localmem(eltype(partials), (256,))
    @inbounds shared[local_i] = i <= count ? abs2(field[i]) : zero(eltype(partials))
    @synchronize
    for stride in (128, 64, 32, 16, 8, 4, 2, 1)
        if local_i <= stride
            @inbounds shared[local_i] += shared[local_i + stride]
        end
        @synchronize
    end
    if local_i == 1
        @inbounds partials[group] = shared[1]
    end
end

KernelAbstractions.@kernel function _compiled_paths_power_finalize_kernel!(power, partials, count)
    local_i = @index(Local, Linear)
    shared = @localmem(eltype(partials), (256,))
    value = zero(eltype(partials))
    # count = cld(length(field), 256), so the final i + 256 is at most
    # count + 256 and remains representable as an Int.
    i = local_i
    while i <= count
        @inbounds value += partials[i]
        i += 256
    end
    @inbounds shared[local_i] = value
    @synchronize
    for stride in (128, 64, 32, 16, 8, 4, 2, 1)
        if local_i <= stride
            @inbounds shared[local_i] += shared[local_i + stride]
        end
        @synchronize
    end
    if local_i == 1
        @inbounds power[1] = shared[1]
    end
end

KernelAbstractions.@kernel function _compiled_paths_normalize_kernel!(field, power, count)
    i = @index(Global, Linear)
    if i <= count
        @inbounds field[i] = iszero(power[1]) ? zero(eltype(field)) : field[i]/sqrt(power[1])
    end
end

function _compiled_paths_normalization_workspace(::Matrix, target, ::Type{T}) where {T}
    return PathNormalizationWorkspace(
        FixedSizeVectorDefault{T}(undef, 0), FixedSizeVectorDefault{T}(undef, 0))
end

function _compiled_paths_normalization_workspace(field, target, ::Type{T}) where {T}
    return PathNormalizationWorkspace(
        _compiled_paths_array(target, T, cld(length(field), 256)),
        _compiled_paths_array(target, T, 1))
end

function _compiled_paths_normalize!(field::Matrix, workspace)
    power=sum(abs2,field)
    if iszero(power)
        fill!(field,zero(eltype(field)))
    else
        field ./= sqrt(power)
    end
    return field
end

function _compiled_paths_normalize!(field, workspace)
    style = Backends.execution_style(field)
    Backends.launch_kernel_async!(style, _compiled_paths_power_partials_kernel!,
        workspace.partials, field, length(field); ndrange=256length(workspace.partials))
    Backends.launch_kernel_async!(style, _compiled_paths_power_finalize_kernel!,
        workspace.power, workspace.partials, length(workspace.partials); ndrange=256)
    Backends.launch_kernel_async!(style, _compiled_paths_normalize_kernel!,
        field, workspace.power, length(field); ndrange=length(field))
    return field
end

# KA masks lanes outside the output axes. Prepared shifts have equal positive
# input/output sizes and 0 <= shift <= axis length. Each shifted coordinate
# lies between 1 and twice the axis length; one subtraction gives a valid index.
KernelAbstractions.@kernel function _compiled_paths_shift_kernel!(out, input, sy, sx)
    y, x = @index(Global, NTuple)
    ny, nx = size(input)
    source_y = y + sy
    source_x = x + sx
    source_y > ny && (source_y -= ny)
    source_x > nx && (source_x -= nx)
    @inbounds out[y, x] = input[source_y, source_x]
end

# Cold preparation bounds crop_start by input_size - output_size; phantom
# crops use nonnegative padding with input_size = output_size + 2padding.
KernelAbstractions.@kernel function _compiled_paths_crop_kernel!(out, input, start)
    y, x = @index(Global, NTuple)
    @inbounds out[y, x] = input[y + start, x + start]
end

# The inside guard proves both source coordinates are within the input axes.
KernelAbstractions.@kernel function _compiled_paths_pad_kernel!(out, input, pad)
    y, x = @index(Global, NTuple)
    inside = pad < y <= size(input, 1) + pad && pad < x <= size(input, 2) + pad
    @inbounds out[y, x] = inside ? input[y - pad, x - pad] : zero(eltype(out))
end

function _compiled_paths_shift!(out, input, inverse::Bool=false)
    n, m = size(input)
    sy = inverse ? fld(n, 2) : cld(n, 2)
    sx = inverse ? fld(m, 2) : cld(m, 2)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_shift_kernel!, out, input, sy, sx; ndrange=size(out))
    return out
end

function _compiled_paths_crop!(out, input, start)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_crop_kernel!, out, input, start; ndrange=size(out))
    return out
end

function _compiled_paths_pad!(out, input, pad)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_pad_kernel!, out, input, pad; ndrange=size(out))
    return out
end


function _compiled_paths_shift!(out::Matrix, input::Matrix, inverse::Bool=false)
    n, m = size(input)
    size(out) == (n, m) || throw(DimensionMismatch("compiled optical paths shift shapes differ"))
    sy = inverse ? fld(n, 2) : cld(n, 2)
    sx = inverse ? fld(m, 2) : cld(m, 2)
    # Both arrays are one-based dense matrices of the checked equal dimensions;
    # each single wrap below maps to a valid input index.
    @inbounds for x in 1:m
        source_x = x + sx
        source_x > m && (source_x -= m)
        for y in 1:n
            source_y = y + sy
            source_y > n && (source_y -= n)
            out[y, x] = input[source_y, source_x]
        end
    end
    return out
end

function _compiled_paths_crop!(out::Matrix, input::Matrix, start)
    for x in axes(out, 2), y in axes(out, 1)
        out[y, x] = input[y + start, x + start]
    end
    return out
end

function _compiled_paths_pad!(out::Matrix, input::Matrix, pad)
    fill!(out, zero(eltype(out)))
    for x in axes(input, 2), y in axes(input, 1)
        out[y + pad, x + pad] = input[y, x]
    end
    return out
end

function _compiled_paths_fft!(out, input, fft_workspace, inverse=false)
    _compiled_paths_shift!(fft_workspace.field, input, true)
    if inverse
        Backends.enqueue_fft_plan!(fft_workspace.field, fft_workspace.inverse)
    else
        Backends.enqueue_fft_plan!(fft_workspace.field, fft_workspace.forward)
    end
    _compiled_paths_shift!(out, fft_workspace.field)
    return out
end

function _compiled_paths_filter!(out, input, transfer, fft_workspace)
    _compiled_paths_fft!(fft_workspace.scratch, input, fft_workspace)
    fft_workspace.scratch .*= transfer
    _compiled_paths_fft!(out, fft_workspace.scratch, fft_workspace, true)
    return out
end

# Pack a centred field into one corner-order FFT plane. Preparation establishes
# size(field,1:2) == size(input) .+ 2pad; the inside test bounds input accesses.
KernelAbstractions.@kernel function _compiled_paths_pack_filter_kernel!(field, input, slot, pad, enabled)
    y, x = @index(Global, NTuple)
    n = size(field, 1)
    sy = y + fld(n, 2)
    sx = x + fld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    inside = pad < sy <= size(input, 1) + pad && pad < sx <= size(input, 2) + pad
    @inbounds field[y, x, slot] = enabled && inside ? input[sy - pad, sx - pad] : zero(eltype(field))
end

# The defocus result remains corner ordered. Read its centred crop while
# packing the next image FFT; this combines crop, fftshift and ifftshift.
KernelAbstractions.@kernel function _compiled_paths_pack_defocus_kernel!(field, defocus, slot, defocus_slot, pad)
    y, x = @index(Global, NTuple)
    n, p = size(field, 1), size(defocus, 1)
    sy = y + fld(n, 2)
    sx = x + fld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    sy += pad + cld(p, 2)
    sx += pad + cld(p, 2)
    sy > p && (sy -= p)
    sx > p && (sx -= p)
    @inbounds field[y, x, slot] = defocus[sy, sx, defocus_slot]
end

# Equal positive square dimensions were checked during preparation. One wrap
# maps each output pixel to the centred inverse-transform result.
KernelAbstractions.@kernel function _compiled_paths_unpack_filter_kernel!(out, field, phase, slot)
    y, x = @index(Global, NTuple)
    n = size(out, 1)
    sy = y + cld(n, 2)
    sx = x + cld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    @inbounds out[y, x] = field[sy, sx, slot] * phase[y, x]
end

function _compiled_paths_pack_filter!(field, input, slot, pad, enabled=true)
    Backends.launch_kernel_async!(Backends.execution_style(field),
        _compiled_paths_pack_filter_kernel!, field, input, slot, pad, enabled;
        ndrange=(size(field, 1), size(field, 2)))
    return nothing
end

function _compiled_paths_pack_defocus!(field, defocus, slot, defocus_slot, pad)
    Backends.launch_kernel_async!(Backends.execution_style(field),
        _compiled_paths_pack_defocus_kernel!, field, defocus, slot, defocus_slot, pad;
        ndrange=(size(field, 1), size(field, 2)))
    return nothing
end

function _compiled_paths_unpack_filter!(out, field, phase, slot)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_unpack_filter_kernel!, out, field, phase, slot; ndrange=size(out))
    return nothing
end

function _compiled_paths_pack_filter!(field::Array, input::Matrix, slot, pad, enabled=true)
    n = size(field, 1)
    for x in 1:n, y in 1:n
        sy, sx = mod1(y + fld(n, 2), n), mod1(x + fld(n, 2), n)
        inside = pad < sy <= size(input, 1) + pad && pad < sx <= size(input, 2) + pad
        field[y, x, slot] = enabled && inside ? input[sy - pad, sx - pad] : zero(eltype(field))
    end
    return nothing
end

function _compiled_paths_pack_defocus!(field::Array, defocus::Array, slot, defocus_slot, pad)
    n, p = size(field, 1), size(defocus, 1)
    for x in 1:n, y in 1:n
        sy = mod1(mod1(y + fld(n, 2), n) + pad + cld(p, 2), p)
        sx = mod1(mod1(x + fld(n, 2), n) + pad + cld(p, 2), p)
        field[y, x, slot] = defocus[sy, sx, defocus_slot]
    end
    return nothing
end

function _compiled_paths_unpack_filter!(out::Matrix, field::Array, phase::Matrix, slot)
    n = size(out, 1)
    for x in 1:n, y in 1:n
        sy, sx = mod1(y + cld(n, 2), n), mod1(x + cld(n, 2), n)
        out[y, x] = field[sy, sx, slot] * phase[y, x]
    end
    return nothing
end

function _compiled_paths_filter_batch!(batch::PathFilterBatch)
    Backends.enqueue_fft_plan!(batch.field, batch.forward)
    batch.field .*= batch.transfer
    Backends.enqueue_fft_plan!(batch.field, batch.inverse)
    return nothing
end
_compiled_paths_filter_batch!(::Nothing) = nothing

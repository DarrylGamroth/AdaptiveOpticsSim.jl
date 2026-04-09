using AdaptiveOpticsSim
using KernelAbstractions

try
    using AMDGPU
catch err
    error("amdgpu_sh_localmem_repro.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

const REPRO_GROUPSIZE = 256

@kernel function sh_cutoff_stats_serial_repro!(stats, spot_cube, valid_mask, cutoff, n_sub::Int, n1::Int, n2::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        base = 3 * (idx - 1)
        T = eltype(stats)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                if val >= cutoff
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
            end
            if total > zero(T)
                sx /= total
                sy /= total
            end
        end
        @inbounds begin
            stats[base + 1] = total
            stats[base + 2] = sx
            stats[base + 3] = sy
        end
    end
end

@kernel function sh_cutoff_stats_localmem_repro!(stats, spot_cube, valid_mask, cutoff, n_sub::Int, n1::Int, n2::Int)
    spot_idx = @index(Group, Linear)
    lid = @index(Local, Linear)
    @uniform group_size = prod(@groupsize())
    n_spots = n_sub * n_sub
    if spot_idx <= n_spots
        i = (spot_idx - 1) ÷ n_sub + 1
        j = spot_idx - (i - 1) * n_sub
        base = 3 * (spot_idx - 1)
        T = eltype(stats)
        totals = @localmem T (group_size,)
        sxs = @localmem T (group_size,)
        sys = @localmem T (group_size,)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            npix = n1 * n2
            p = lid
            while p <= npix
                x = fld(p - 1, n2) + 1
                y = mod(p - 1, n2) + 1
                val = @inbounds spot_cube[spot_idx, x, y]
                if val >= cutoff
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
                p += group_size
            end
        end
        @inbounds begin
            totals[lid] = total
            sxs[lid] = sx
            sys[lid] = sy
        end
        @synchronize

        stride = group_size >>> 1
        while stride > 0
            if lid <= stride
                @inbounds begin
                    totals[lid] += totals[lid + stride]
                    sxs[lid] += sxs[lid + stride]
                    sys[lid] += sys[lid + stride]
                end
            end
            @synchronize
            stride >>>= 1
        end

        if lid == 1
            total_sum = @inbounds totals[1]
            sx_sum = @inbounds sxs[1]
            sy_sum = @inbounds sys[1]
            if total_sum > zero(T)
                sx_sum /= total_sum
                sy_sum /= total_sum
            end
            @inbounds begin
                stats[base + 1] = total_sum
                stats[base + 2] = sx_sum
                stats[base + 3] = sy_sum
            end
        end
    end
end

function build_repro_inputs(::Type{T}=Float32) where {T<:AbstractFloat}
    n_sub = 2
    n1 = 22
    n2 = 22
    n_spots = n_sub * n_sub
    spot_cube = reshape(T.(1:(n_spots * n1 * n2)), n_spots, n1, n2)
    valid_mask = trues(n_sub, n_sub)
    stats = zeros(T, 3 * n_spots)
    cutoff = T(0.1)
    return (; stats, spot_cube, valid_mask, cutoff, n_sub, n1, n2, n_spots)
end

function run_serial_repro()
    data = build_repro_inputs()
    backend = AMDGPU.ROCBackend()
    stats = AMDGPU.ROCArray(data.stats)
    spot_cube = AMDGPU.ROCArray(data.spot_cube)
    valid_mask = AMDGPU.ROCArray(data.valid_mask)
    kernel! = sh_cutoff_stats_serial_repro!(backend)
    kernel!(stats, spot_cube, valid_mask, data.cutoff, data.n_sub, data.n1, data.n2; ndrange=data.n_spots)
    KernelAbstractions.synchronize(backend)
    return Array(stats)
end

function run_localmem_repro()
    data = build_repro_inputs()
    backend = AMDGPU.ROCBackend()
    stats = AMDGPU.ROCArray(data.stats)
    spot_cube = AMDGPU.ROCArray(data.spot_cube)
    valid_mask = AMDGPU.ROCArray(data.valid_mask)
    kernel! = sh_cutoff_stats_localmem_repro!(backend, REPRO_GROUPSIZE)
    kernel!(stats, spot_cube, valid_mask, data.cutoff, data.n_sub, data.n1, data.n2; ndrange=data.n_spots * REPRO_GROUPSIZE)
    KernelAbstractions.synchronize(backend)
    return Array(stats)
end

mode = length(ARGS) == 0 ? "localmem" : ARGS[1]
if mode == "serial"
    stats = run_serial_repro()
    println((mode=mode, sum=sum(stats), first=stats[1:6]))
elseif mode == "localmem"
    stats = run_localmem_repro()
    println((mode=mode, sum=sum(stats), first=stats[1:6]))
else
    error("unsupported mode '$mode'; use serial or localmem")
end

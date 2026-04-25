using AdaptiveOpticsSim
using Random
using DelimitedFiles
using SparseArrays
using LinearAlgebra: mul!
using KernelAbstractions: @kernel, @index

@kernel function revolt_scatter_active_command_kernel!(full_command, active_command, active_indices, n_active::Int)
    i = @index(Global, Linear)
    if i <= n_active
        @inbounds full_command[active_indices[i]] = active_command[i]
    end
end

@kernel function revolt_sh_spot_mosaic_kernel!(mosaic, spot_cube, n_lenslets::Int, roi::Int)
    idx, u, v = @index(Global, NTuple)
    if idx <= size(spot_cube, 1) && u <= roi && v <= roi
        sub_i = mod1(idx, n_lenslets)
        sub_j = fld(idx - 1, n_lenslets) + 1
        out_i = (sub_i - 1) * roi + u
        out_j = (sub_j - 1) * roi + v
        @inbounds mosaic[out_i, out_j] = spot_cube[idx, u, v]
    end
end

struct RevoltLikeHILContext{TEL,SRC,DM,WFS,DET,AI,AC,EC,EB,TF,RNG}
    tel::TEL
    src::SRC
    dm::DM
    wfs::WFS
    det::DET
    active_indices_backend::AI
    active_command::AC
    extrapolated_command::EC
    extrapolation_backend::EB
    tiled_frame::TF
    n_lenslets::Int
    roi::Int
    rng::RNG
end

function revolt_profile_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return (; selector=CPUBackend(), array_backend=Array, label="cpu")
    elseif lowered == "cuda"
        @eval import CUDA
        CUDA.functional() || error("REVOLT-like HIL profiling requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        array_backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        array_backend === nothing && error("CUDA backend array type is unavailable")
        return (; selector=CUDABackend(), array_backend, label="cuda")
    elseif lowered == "amdgpu"
        @eval import AMDGPU
        AMDGPU.functional() || error("REVOLT-like HIL profiling requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        array_backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        array_backend === nothing && error("AMDGPU backend array type is unavailable")
        return (; selector=AMDGPUBackend(), array_backend, label="amdgpu")
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function revolt_load_dm277_actuator_map(path::AbstractString)
    raw = readlines(path)
    length(raw) >= 2 || error("actuator map file is unexpectedly short")
    header = split(strip(raw[1]))
    length(header) == 4 || error("unexpected actuator map header '$((raw[1]))'")
    n_row = parse(Int, header[3])
    n_col = parse(Int, header[4])
    grid = DelimitedFiles.readdlm(IOBuffer(join(raw[2:end], '\n')), ',', Int)
    size(grid) == (n_row, n_col) || error("actuator map grid size $(size(grid)) does not match header $((n_row, n_col))")
    maximum(grid) > 0 || error("actuator map has no active actuators")
    active_indices = Vector{Int}(undef, maximum(grid))
    @inbounds for j in 1:n_col, i in 1:n_row
        label = grid[i, j]
        if label > 0
            active_indices[label] = LinearIndices(grid)[i, j]
        end
    end
    return grid, active_indices
end

function revolt_load_dm_extrapolation(path::AbstractString, ::Type{T}) where {T<:AbstractFloat}
    raw = readlines(path)
    length(raw) >= 2 || error("extrapolation file is unexpectedly short")
    header = split(strip(raw[1]))
    length(header) == 4 || error("unexpected extrapolation header '$((raw[1]))'")
    n_row = parse(Int, header[3])
    n_col = parse(Int, header[4])
    triplets = DelimitedFiles.readdlm(IOBuffer(join(raw[2:end], '\n')), ',', Float64)
    size(triplets, 2) == 3 || error("dmExtrapolation triplets must have three columns")
    rows = Int.(triplets[:, 1]) .+ 1
    cols = Int.(triplets[:, 2]) .+ 1
    vals = T.(triplets[:, 3])
    return SparseArrays.sparse(rows, cols, vals, n_row, n_col)
end

function revolt_fill_active_command!(active_command::AbstractVector{T}) where {T<:AbstractFloat}
    host = Vector{T}(undef, length(active_command))
    @inbounds for i in eachindex(host)
        host[i] = T(sin(0.013 * i) + 0.25 * cos(0.031 * i))
    end
    copyto!(active_command, host)
    return active_command
end

function revolt_scatter_active_command!(full_command::AbstractVector{T},
    active_command::AbstractVector{T}, active_indices_backend::AbstractVector{Int}) where {T<:AbstractFloat}
    revolt_scatter_active_command!(AdaptiveOpticsSim.execution_style(full_command), full_command, active_command, active_indices_backend)
    return full_command
end

function revolt_scatter_active_command!(::AdaptiveOpticsSim.ScalarCPUStyle, full_command::AbstractVector{T},
    active_command::AbstractVector{T}, active_indices_backend::AbstractVector{Int}) where {T<:AbstractFloat}
    fill!(full_command, zero(T))
    @inbounds for i in eachindex(active_command)
        full_command[active_indices_backend[i]] = active_command[i]
    end
    return full_command
end

function revolt_scatter_active_command!(style::AdaptiveOpticsSim.AcceleratorStyle, full_command::AbstractVector{T},
    active_command::AbstractVector{T}, active_indices_backend::AbstractVector{Int}) where {T<:AbstractFloat}
    fill!(full_command, zero(T))
    AdaptiveOpticsSim.launch_kernel!(style, revolt_scatter_active_command_kernel!,
        full_command, active_command, active_indices_backend, length(active_command);
        ndrange=length(active_command))
    return full_command
end

function revolt_tile_spot_cube!(mosaic::AbstractMatrix{T}, spot_cube::AbstractArray{T,3}, n_lenslets::Int, roi::Int) where {T<:AbstractFloat}
    revolt_tile_spot_cube!(AdaptiveOpticsSim.execution_style(mosaic), mosaic, spot_cube, n_lenslets, roi)
    return mosaic
end

function revolt_tile_spot_cube!(::AdaptiveOpticsSim.ScalarCPUStyle, mosaic::AbstractMatrix{T},
    spot_cube::AbstractArray{T,3}, n_lenslets::Int, roi::Int) where {T<:AbstractFloat}
    @inbounds for idx in 1:size(spot_cube, 1)
        sub_i = mod1(idx, n_lenslets)
        sub_j = fld(idx - 1, n_lenslets) + 1
        out_i = (sub_i - 1) * roi + 1
        out_j = (sub_j - 1) * roi + 1
        copyto!(@view(mosaic[out_i:out_i + roi - 1, out_j:out_j + roi - 1]), @view(spot_cube[idx, :, :]))
    end
    return mosaic
end

function revolt_tile_spot_cube!(style::AdaptiveOpticsSim.AcceleratorStyle, mosaic::AbstractMatrix{T},
    spot_cube::AbstractArray{T,3}, n_lenslets::Int, roi::Int) where {T<:AbstractFloat}
    AdaptiveOpticsSim.launch_kernel!(style, revolt_sh_spot_mosaic_kernel!, mosaic, spot_cube, n_lenslets, roi;
        ndrange=size(spot_cube))
    return mosaic
end

function build_revolt_like_hil_context(; backend_name::AbstractString="cpu", config_dir::AbstractString,
    sensor=CMOSSensor(), response_model=nothing,
    thermal_model=nothing, T::Type{<:AbstractFloat}=Float32,
    rng=MersenneTwister(0))
    backend_cfg = revolt_profile_backend(backend_name)
    dark_current = isnothing(thermal_model) ? zero(T) : T(0.02)
    actuator_map_path = joinpath(config_dir, "revolt_like_dmActuatorMap_277.csv")
    extrapolation_path = joinpath(config_dir, "revolt_like_dmExtrapolation.csv")
    actuator_map, active_indices_host = revolt_load_dm277_actuator_map(actuator_map_path)
    extrapolation_host = Matrix(revolt_load_dm_extrapolation(extrapolation_path, T))
    n_lenslets = 16
    roi = 22
    resolution = 352
    n_act = size(actuator_map, 1)
    n_active = length(active_indices_host)

    tel = Telescope(
        resolution=resolution,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        T=T,
        backend=backend_cfg.selector,
    )
    src = Source(band=:I, magnitude=0.0, T=T)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.3, T=T, backend=backend_cfg.selector)
    wfs = ShackHartmann(tel; n_lenslets=n_lenslets, mode=Diffractive(), n_pix_subap=roi,
        diffraction_padding=2, T=T, backend=backend_cfg.selector)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1,
        dark_current=dark_current, sensor=sensor, response_model=response_model,
        thermal_model=thermal_model, T=T, backend=backend_cfg.selector)
    active_indices_backend = backend_cfg.array_backend{Int}(undef, n_active)
    copyto!(active_indices_backend, active_indices_host)
    active_command = backend_cfg.array_backend{T}(undef, n_active)
    extrapolated_command = backend_cfg.array_backend{T}(undef, n_active)
    extrapolation_backend = backend_cfg.array_backend{T}(undef, size(extrapolation_host)...)
    copyto!(extrapolation_backend, extrapolation_host)
    revolt_fill_active_command!(active_command)
    tiled_frame = backend_cfg.array_backend{T}(undef, resolution, resolution)

    AdaptiveOpticsSim.ensure_sh_calibration!(wfs, tel, src)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs.state.slopes))

    return RevoltLikeHILContext(
        tel,
        src,
        dm,
        wfs,
        det,
        active_indices_backend,
        active_command,
        extrapolated_command,
        extrapolation_backend,
        tiled_frame,
        n_lenslets,
        roi,
        rng,
    )
end

function revolt_like_command_map!(ctx::RevoltLikeHILContext)
    reset_opd!(ctx.tel)
    mul!(ctx.extrapolated_command, ctx.extrapolation_backend, ctx.active_command)
    revolt_scatter_active_command!(ctx.dm.state.coefs, ctx.extrapolated_command, ctx.active_indices_backend)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.dm.state.coefs))
    return nothing
end

function revolt_like_dm_apply!(ctx::RevoltLikeHILContext)
    mul!(ctx.extrapolated_command, ctx.extrapolation_backend, ctx.active_command)
    revolt_scatter_active_command!(ctx.dm.state.coefs, ctx.extrapolated_command, ctx.active_indices_backend)
    apply!(ctx.dm, ctx.tel, DMReplace())
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.tel.state.opd))
    return nothing
end

function revolt_like_sense!(ctx::RevoltLikeHILContext)
    mul!(ctx.extrapolated_command, ctx.extrapolation_backend, ctx.active_command)
    revolt_scatter_active_command!(ctx.dm.state.coefs, ctx.extrapolated_command, ctx.active_indices_backend)
    apply!(ctx.dm, ctx.tel, DMReplace())
    measure!(ctx.wfs, ctx.tel, ctx.src, ctx.det; rng=ctx.rng)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.wfs.state.spot_cube))
    return nothing
end

function revolt_like_mosaic!(ctx::RevoltLikeHILContext)
    mul!(ctx.extrapolated_command, ctx.extrapolation_backend, ctx.active_command)
    revolt_scatter_active_command!(ctx.dm.state.coefs, ctx.extrapolated_command, ctx.active_indices_backend)
    apply!(ctx.dm, ctx.tel, DMReplace())
    measure!(ctx.wfs, ctx.tel, ctx.src, ctx.det; rng=ctx.rng)
    revolt_tile_spot_cube!(ctx.tiled_frame, ctx.wfs.state.spot_cube, ctx.n_lenslets, ctx.roi)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.tiled_frame))
    return nothing
end

function revolt_like_step!(ctx::RevoltLikeHILContext)
    reset_opd!(ctx.tel)
    mul!(ctx.extrapolated_command, ctx.extrapolation_backend, ctx.active_command)
    revolt_scatter_active_command!(ctx.dm.state.coefs, ctx.extrapolated_command, ctx.active_indices_backend)
    apply!(ctx.dm, ctx.tel, DMReplace())
    measure!(ctx.wfs, ctx.tel, ctx.src, ctx.det; rng=ctx.rng)
    revolt_tile_spot_cube!(ctx.tiled_frame, ctx.wfs.state.spot_cube, ctx.n_lenslets, ctx.roi)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.dm.state.coefs))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.wfs.state.spot_cube))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(ctx.tiled_frame))
    return nothing
end

function revolt_like_allocation_profile(; backend_name::AbstractString="cpu", config_dir::AbstractString,
    sensor=CMOSSensor(), response_model=nothing,
    thermal_model=nothing, T::Type{<:AbstractFloat}=Float32)
    ctx = build_revolt_like_hil_context(; backend_name, config_dir, sensor, response_model, thermal_model, T)
    revolt_like_command_map!(ctx)
    revolt_like_dm_apply!(ctx)
    revolt_like_sense!(ctx)
    revolt_like_mosaic!(ctx)
    revolt_like_step!(ctx)
    GC.gc()
    return (;
        command_map=(@allocated revolt_like_command_map!(ctx)),
        dm_apply=(@allocated revolt_like_dm_apply!(ctx)),
        sense=(@allocated revolt_like_sense!(ctx)),
        mosaic=(@allocated revolt_like_mosaic!(ctx)),
        total=(@allocated revolt_like_step!(ctx)),
    )
end

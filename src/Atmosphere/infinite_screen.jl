@kernel function gather_stencil_data_kernel!(dest, screen, coords, n::Int)
    k = @index(Global, Linear)
    if k <= n
        @inbounds dest[k] = screen[coords[k, 2], coords[k, 1]]
    end
end

@kernel function inject_column_positive_kernel!(dest, src, boundary, m::Int)
    i, j = @index(Global, NTuple)
    if i <= m && j <= m
        @inbounds dest[i, j] = j == 1 ? boundary[i] : src[i, j - 1]
    end
end

@kernel function inject_column_negative_kernel!(dest, src, boundary, m::Int)
    i, j = @index(Global, NTuple)
    if i <= m && j <= m
        @inbounds dest[i, j] = j == m ? boundary[i] : src[i, j + 1]
    end
end

@kernel function inject_row_positive_kernel!(dest, src, boundary, m::Int)
    i, j = @index(Global, NTuple)
    if i <= m && j <= m
        @inbounds dest[i, j] = i == 1 ? boundary[j] : src[i - 1, j]
    end
end

@kernel function inject_row_negative_kernel!(dest, src, boundary, m::Int)
    i, j = @index(Global, NTuple)
    if i <= m && j <= m
        @inbounds dest[i, j] = i == m ? boundary[j] : src[i + 1, j]
    end
end

"""
Builder-time parameters for a single infinite phase screen.

The maintained screen extracted under the pupil has resolution
`screen_resolution`, while the persistent internal screen buffer and injected
boundary length both use `stencil_size`.
"""
struct InfinitePhaseScreenParams{T<:AbstractFloat}
    r0::T
    L0::T
    pixel_scale::T
    screen_resolution::Int
    stencil_size::Int
end

"""
Mutable runtime state for a single infinite phase screen.

- `screen` is the persistent full-screen buffer of size `stencil_size`.
- `extract_buffer` is the pupil-sized sampling workspace.
- `column_positive` / `column_negative` / `row_positive` / `row_negative`
  hold the precomputed boundary-injection models for each screen edge.
"""
mutable struct InfinitePhaseScreenState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    V<:AbstractVector{T},
    I<:AbstractMatrix{Int},
    B,
}
    screen::A
    screen_scratch::A
    extract_buffer::A
    stencil_buffer::V
    boundary_buffer::V
    noise_buffer::V
    column_positive_coords::I
    column_negative_coords::I
    row_positive_coords::I
    row_negative_coords::I
    column_positive::B
    column_negative::B
    row_positive::B
    row_negative::B
    initialized::Bool
end

"""
Container for a single infinite phase screen and its state.
"""
struct InfinitePhaseScreen{
    P<:InfinitePhaseScreenParams,
    S<:InfinitePhaseScreenState,
    G<:KolmogorovAtmosphere,
    TT<:Telescope,
}
    params::P
    state::S
    generator::G
    generator_telescope::TT
end

"""
Per-layer metadata for the infinite multilayer atmosphere backend.
"""
struct InfiniteLayerParams{T<:AbstractFloat}
    amplitude_scale::T
    wind_velocity_x::T
    wind_velocity_y::T
    altitude::T
end

"""
Per-layer transport state for the infinite multilayer atmosphere backend.
"""
mutable struct InfiniteLayerState{T<:AbstractFloat}
    offset_x::T
    offset_y::T
    integer_shift_x::Int
    integer_shift_y::Int
end

"""
One atmospheric layer in the planned infinite boundary-injection model.

The future runtime will reuse the existing pupil extraction helper from
`multilayer.jl` after updating the persistent screen buffer through boundary
injection instead of periodic wraparound.
"""
struct InfiniteAtmosphereLayer{
    P<:InfiniteLayerParams,
    S<:InfiniteLayerState,
    Screen<:InfinitePhaseScreen,
}
    params::P
    screen::Screen
    state::S
end

"""
Parameters for the planned infinite multilayer atmosphere backend.
"""
struct InfiniteMultiLayerParams{T<:AbstractFloat,
    V1<:AbstractVector{T},
    V2<:AbstractVector{T},
    V3<:AbstractVector{T},
    V4<:AbstractVector{T},
    V5<:AbstractVector{T},
    V6<:AbstractVector{T}}
    cn2_fractions::V1
    wind_speed::V2
    wind_direction::V3
    altitude::V4
    wind_velocity_x::V5
    wind_velocity_y::V6
    r0::T
    L0::T
end

"""
Runtime state for the planned infinite multilayer atmosphere backend.
"""
mutable struct InfiniteMultiLayerState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
    layer_buffer::A
end

"""
Infinite boundary-injection atmosphere backend.

This backend supports both CPU and GPU runtime stepping. Builder-side
factorization remains separate from the steady-state transport path so the hot
loop stays backend-generic.
"""
struct InfiniteMultiLayerAtmosphere{
    P<:InfiniteMultiLayerParams,
    S<:InfiniteMultiLayerState,
    L<:AbstractVector{<:InfiniteAtmosphereLayer},
} <: AbstractAtmosphere
    params::P
    layers::L
    state::S
end

@inline default_infinite_screen_resolution(n::Int) = 3 * n

@inline default_infinite_stencil_size(n::Int) = max(257, 2 * default_infinite_screen_resolution(n) + 1)

function infinite_screen_telescope(
    tel::Telescope;
    resolution::Int,
    T::Type{<:AbstractFloat}=Float64,
    backend=Array,
)
    delta = tel.params.diameter / tel.params.resolution
    return Telescope(
        resolution=resolution,
        diameter=delta * resolution,
        sampling_time=tel.params.sampling_time,
        central_obstruction=0.0,
        fov_arcsec=tel.params.fov_arcsec,
        T=T,
        backend=backend,
    )
end

function _fill_stencil_data!(dest::AbstractVector{T}, screen::AbstractMatrix{T},
    model::InfiniteBoundaryModel) where {T<:AbstractFloat}
    coords = model.stencil.stencil_coords
    length(dest) == size(coords, 1) ||
        throw(DimensionMismatchError("stencil buffer length must match the boundary-model stencil length"))
    @inbounds for k in eachindex(dest)
        dest[k] = screen[coords[k, 2], coords[k, 1]]
    end
    return dest
end

function _fill_stencil_data!(style::AcceleratorStyle, dest::AbstractVector{T}, screen::AbstractMatrix{T},
    coords::AbstractMatrix{Int}) where {T<:AbstractFloat}
    launch_kernel!(style, gather_stencil_data_kernel!, dest, screen, coords, length(dest); ndrange=length(dest))
    return dest
end

@inline function _swap_screen_buffers!(state::InfinitePhaseScreenState)
    screen = state.screen
    state.screen = state.screen_scratch
    state.screen_scratch = screen
    return state
end

function _inject_column_positive!(screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    model = state.column_positive
    _fill_stencil_data!(state.stencil_buffer, screen, model)
    sample_boundary_line!(state.boundary_buffer, model.operator, state.stencil_buffer, state.noise_buffer, rng)
    @inbounds for row in 1:m
        for col in m:-1:2
            screen[row, col] = screen[row, col - 1]
        end
        screen[row, 1] = state.boundary_buffer[row]
    end
    return screen
end

@inline _inject_column_positive!(::ScalarCPUStyle, screen::AbstractMatrix{T},
    state::InfinitePhaseScreenState{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    _inject_column_positive!(screen, state, rng)

function _inject_column_positive!(style::AcceleratorStyle, screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    _fill_stencil_data!(style, state.stencil_buffer, screen, state.column_positive_coords)
    sample_boundary_line!(state.boundary_buffer, state.column_positive.operator, state.stencil_buffer, state.noise_buffer, rng)
    launch_kernel!(style, inject_column_positive_kernel!, state.screen_scratch, screen, state.boundary_buffer, m; ndrange=size(screen))
    _swap_screen_buffers!(state)
    return state.screen
end

function _inject_column_negative!(screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    model = state.column_negative
    _fill_stencil_data!(state.stencil_buffer, screen, model)
    sample_boundary_line!(state.boundary_buffer, model.operator, state.stencil_buffer, state.noise_buffer, rng)
    @inbounds for row in 1:m
        for col in 1:(m - 1)
            screen[row, col] = screen[row, col + 1]
        end
        screen[row, m] = state.boundary_buffer[row]
    end
    return screen
end

@inline _inject_column_negative!(::ScalarCPUStyle, screen::AbstractMatrix{T},
    state::InfinitePhaseScreenState{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    _inject_column_negative!(screen, state, rng)

function _inject_column_negative!(style::AcceleratorStyle, screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    _fill_stencil_data!(style, state.stencil_buffer, screen, state.column_negative_coords)
    sample_boundary_line!(state.boundary_buffer, state.column_negative.operator, state.stencil_buffer, state.noise_buffer, rng)
    launch_kernel!(style, inject_column_negative_kernel!, state.screen_scratch, screen, state.boundary_buffer, m; ndrange=size(screen))
    _swap_screen_buffers!(state)
    return state.screen
end

function _inject_row_positive!(screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    model = state.row_positive
    _fill_stencil_data!(state.stencil_buffer, screen, model)
    sample_boundary_line!(state.boundary_buffer, model.operator, state.stencil_buffer, state.noise_buffer, rng)
    @inbounds begin
        for row in m:-1:2
            for col in 1:m
                screen[row, col] = screen[row - 1, col]
            end
        end
        for col in 1:m
            screen[1, col] = state.boundary_buffer[col]
        end
    end
    return screen
end

@inline _inject_row_positive!(::ScalarCPUStyle, screen::AbstractMatrix{T},
    state::InfinitePhaseScreenState{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    _inject_row_positive!(screen, state, rng)

function _inject_row_positive!(style::AcceleratorStyle, screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    _fill_stencil_data!(style, state.stencil_buffer, screen, state.row_positive_coords)
    sample_boundary_line!(state.boundary_buffer, state.row_positive.operator, state.stencil_buffer, state.noise_buffer, rng)
    launch_kernel!(style, inject_row_positive_kernel!, state.screen_scratch, screen, state.boundary_buffer, m; ndrange=size(screen))
    _swap_screen_buffers!(state)
    return state.screen
end

function _inject_row_negative!(screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    model = state.row_negative
    _fill_stencil_data!(state.stencil_buffer, screen, model)
    sample_boundary_line!(state.boundary_buffer, model.operator, state.stencil_buffer, state.noise_buffer, rng)
    @inbounds begin
        for row in 1:(m - 1)
            for col in 1:m
                screen[row, col] = screen[row + 1, col]
            end
        end
        for col in 1:m
            screen[m, col] = state.boundary_buffer[col]
        end
    end
    return screen
end

@inline _inject_row_negative!(::ScalarCPUStyle, screen::AbstractMatrix{T},
    state::InfinitePhaseScreenState{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    _inject_row_negative!(screen, state, rng)

function _inject_row_negative!(style::AcceleratorStyle, screen::AbstractMatrix{T}, state::InfinitePhaseScreenState{T},
    rng::AbstractRNG) where {T<:AbstractFloat}
    m = size(screen, 1)
    _fill_stencil_data!(style, state.stencil_buffer, screen, state.row_negative_coords)
    sample_boundary_line!(state.boundary_buffer, state.row_negative.operator, state.stencil_buffer, state.noise_buffer, rng)
    launch_kernel!(style, inject_row_negative_kernel!, state.screen_scratch, screen, state.boundary_buffer, m; ndrange=size(screen))
    _swap_screen_buffers!(state)
    return state.screen
end

function ensure_initialized!(screen::InfinitePhaseScreen, rng::AbstractRNG)
    if !screen.state.initialized
        advance!(screen.generator, screen.generator_telescope, rng)
        copyto!(screen.state.screen, screen.generator.state.opd)
        screen.state.initialized = true
    end
    return screen
end

function _apply_integer_shift_x!(screen::InfinitePhaseScreen, shift::Int, rng::AbstractRNG)
    style = execution_style(screen.state.screen)
    if shift > 0
        for _ in 1:shift
            _inject_column_positive!(style, screen.state.screen, screen.state, rng)
        end
    elseif shift < 0
        for _ in 1:(-shift)
            _inject_column_negative!(style, screen.state.screen, screen.state, rng)
        end
    end
    return screen
end

function _apply_integer_shift_y!(screen::InfinitePhaseScreen, shift::Int, rng::AbstractRNG)
    style = execution_style(screen.state.screen)
    if shift > 0
        for _ in 1:shift
            _inject_row_positive!(style, screen.state.screen, screen.state, rng)
        end
    elseif shift < 0
        for _ in 1:(-shift)
            _inject_row_negative!(style, screen.state.screen, screen.state, rng)
        end
    end
    return screen
end

function sample_layer!(out::AbstractMatrix{T}, layer::InfiniteAtmosphereLayer, tel::Telescope,
    rng::AbstractRNG) where {T<:AbstractFloat}
    ensure_initialized!(layer.screen, rng)
    delta = T(tel.params.diameter / tel.params.resolution)
    dt = T(tel.params.sampling_time)
    next_offset_x = layer.state.offset_x + T(layer.params.wind_velocity_x) * dt / delta
    next_offset_y = layer.state.offset_y + T(layer.params.wind_velocity_y) * dt / delta
    shift_x = trunc(Int, next_offset_x)
    shift_y = trunc(Int, next_offset_y)
    residual_x = next_offset_x - T(shift_x)
    residual_y = next_offset_y - T(shift_y)
    _apply_integer_shift_x!(layer.screen, shift_x, rng)
    _apply_integer_shift_y!(layer.screen, shift_y, rng)
    layer.state.offset_x = residual_x
    layer.state.offset_y = residual_y
    layer.state.integer_shift_x += shift_x
    layer.state.integer_shift_y += shift_y
    return render_layer!(out, layer, tel)
end

function render_layer!(out::AbstractMatrix{T}, layer::InfiniteAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    shift_x, shift_y, footprint_scale = layer_source_geometry(src, layer.params.altitude, tel, T)
    extract_shifted_screen!(out, layer.screen.state.screen,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        T(layer.params.amplitude_scale),
        footprint_scale)
    return out
end

function _warmup_gpu_infinite_screen!(style::AcceleratorStyle, screen::InfinitePhaseScreen, rng::AbstractRNG)
    state = screen.state
    fill!(state.screen, zero(eltype(state.screen)))
    fill!(state.screen_scratch, zero(eltype(state.screen_scratch)))
    fill!(state.extract_buffer, zero(eltype(state.extract_buffer)))
    fill!(state.stencil_buffer, zero(eltype(state.stencil_buffer)))
    fill!(state.boundary_buffer, zero(eltype(state.boundary_buffer)))
    fill!(state.noise_buffer, zero(eltype(state.noise_buffer)))
    _inject_column_positive!(style, state.screen, state, rng)
    _inject_column_negative!(style, state.screen, state, rng)
    _inject_row_positive!(style, state.screen, state, rng)
    _inject_row_negative!(style, state.screen, state, rng)
    fill!(state.screen, zero(eltype(state.screen)))
    fill!(state.screen_scratch, zero(eltype(state.screen_scratch)))
    fill!(state.extract_buffer, zero(eltype(state.extract_buffer)))
    fill!(state.stencil_buffer, zero(eltype(state.stencil_buffer)))
    fill!(state.boundary_buffer, zero(eltype(state.boundary_buffer)))
    fill!(state.noise_buffer, zero(eltype(state.noise_buffer)))
    synchronize_backend!(style)
    return screen
end

@inline _warmup_gpu_infinite_screen!(::ScalarCPUStyle, screen::InfinitePhaseScreen, rng::AbstractRNG) = screen

function InfinitePhaseScreen(tel::Telescope;
    r0::Real,
    L0::Real=25.0,
    screen_resolution::Int=default_infinite_screen_resolution(tel.params.resolution),
    stencil_size::Int=default_infinite_stencil_size(tel.params.resolution),
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)
    screen_resolution >= tel.params.resolution ||
        throw(InvalidConfiguration("screen_resolution must be >= telescope resolution"))
    pixel_scale = T(tel.params.diameter / tel.params.resolution)
    validate_infinite_stencil(stencil_size, screen_resolution, :column, :positive, screen_resolution)
    params = InfinitePhaseScreenParams(T(r0), T(L0), pixel_scale, screen_resolution, stencil_size)
    column_positive = boundary_injection_model(screen_resolution, pixel_scale, params.r0, params.L0;
        stencil_size=stencil_size,
        orientation=:column,
        side=:positive,
        tail_stride=screen_resolution,
        T=T)
    column_negative = boundary_injection_model(screen_resolution, pixel_scale, params.r0, params.L0;
        stencil_size=stencil_size,
        orientation=:column,
        side=:negative,
        tail_stride=screen_resolution,
        T=T)
    row_positive = boundary_injection_model(screen_resolution, pixel_scale, params.r0, params.L0;
        stencil_size=stencil_size,
        orientation=:row,
        side=:positive,
        tail_stride=screen_resolution,
        T=T)
    row_negative = boundary_injection_model(screen_resolution, pixel_scale, params.r0, params.L0;
        stencil_size=stencil_size,
        orientation=:row,
        side=:negative,
        tail_stride=screen_resolution,
        T=T)
    screen_telescope = infinite_screen_telescope(tel; resolution=stencil_size, T=T, backend=backend)
    generator = KolmogorovAtmosphere(screen_telescope; r0=params.r0, L0=params.L0, T=T, backend=backend)
    screen = backend{T}(undef, stencil_size, stencil_size)
    screen_scratch = backend{T}(undef, stencil_size, stencil_size)
    extract_buffer = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    build_backend = default_build_backend(screen)
    fill!(screen, zero(T))
    fill!(screen_scratch, zero(T))
    fill!(extract_buffer, zero(T))
    state = InfinitePhaseScreenState(
        screen,
        screen_scratch,
        extract_buffer,
        materialize_build(build_backend, backend{T}(undef, size(column_positive.operator.predictor, 2)),
            Vector{T}(undef, size(column_positive.operator.predictor, 2))),
        materialize_build(build_backend, backend{T}(undef, size(column_positive.operator.predictor, 1)),
            Vector{T}(undef, size(column_positive.operator.predictor, 1))),
        materialize_build(build_backend, backend{T}(undef, size(column_positive.operator.residual_factor, 2)),
            Vector{T}(undef, size(column_positive.operator.residual_factor, 2))),
        materialize_build(build_backend, Matrix{Int}(undef, size(column_positive.stencil.stencil_coords)...),
            column_positive.stencil.stencil_coords),
        materialize_build(build_backend, Matrix{Int}(undef, size(column_negative.stencil.stencil_coords)...),
            column_negative.stencil.stencil_coords),
        materialize_build(build_backend, Matrix{Int}(undef, size(row_positive.stencil.stencil_coords)...),
            row_positive.stencil.stencil_coords),
        materialize_build(build_backend, Matrix{Int}(undef, size(row_negative.stencil.stencil_coords)...),
            row_negative.stencil.stencil_coords),
        InfiniteBoundaryModel(
            InfiniteBoundaryStencil(
                column_positive.stencil.stencil_coords,
                column_positive.stencil.boundary_coords,
                column_positive.stencil.stencil_positions,
                column_positive.stencil.boundary_positions,
                column_positive.stencil.orientation,
                column_positive.stencil.side,
            ),
            InfiniteBoundaryOperator(
                materialize_build(build_backend, screen, column_positive.operator.predictor),
                materialize_build(build_backend, screen, column_positive.operator.residual_factor),
                materialize_build(build_backend, screen, column_positive.operator.cov_zz),
                materialize_build(build_backend, screen, column_positive.operator.cov_xx),
                materialize_build(build_backend, screen, column_positive.operator.cov_xz),
                materialize_build(build_backend, screen, column_positive.operator.cov_zx),
                materialize_build(build_backend, screen, column_positive.operator.residual_covariance),
                materialize_build(build_backend, Vector{T}(undef, length(column_positive.operator.singular_values)),
                    column_positive.operator.singular_values),
                T(column_positive.operator.condition_ratio),
                column_positive.operator.orientation,
                column_positive.operator.side,
            ),
        ),
        InfiniteBoundaryModel(
            InfiniteBoundaryStencil(
                column_negative.stencil.stencil_coords,
                column_negative.stencil.boundary_coords,
                column_negative.stencil.stencil_positions,
                column_negative.stencil.boundary_positions,
                column_negative.stencil.orientation,
                column_negative.stencil.side,
            ),
            InfiniteBoundaryOperator(
                materialize_build(build_backend, screen, column_negative.operator.predictor),
                materialize_build(build_backend, screen, column_negative.operator.residual_factor),
                materialize_build(build_backend, screen, column_negative.operator.cov_zz),
                materialize_build(build_backend, screen, column_negative.operator.cov_xx),
                materialize_build(build_backend, screen, column_negative.operator.cov_xz),
                materialize_build(build_backend, screen, column_negative.operator.cov_zx),
                materialize_build(build_backend, screen, column_negative.operator.residual_covariance),
                materialize_build(build_backend, Vector{T}(undef, length(column_negative.operator.singular_values)),
                    column_negative.operator.singular_values),
                T(column_negative.operator.condition_ratio),
                column_negative.operator.orientation,
                column_negative.operator.side,
            ),
        ),
        InfiniteBoundaryModel(
            InfiniteBoundaryStencil(
                row_positive.stencil.stencil_coords,
                row_positive.stencil.boundary_coords,
                row_positive.stencil.stencil_positions,
                row_positive.stencil.boundary_positions,
                row_positive.stencil.orientation,
                row_positive.stencil.side,
            ),
            InfiniteBoundaryOperator(
                materialize_build(build_backend, screen, row_positive.operator.predictor),
                materialize_build(build_backend, screen, row_positive.operator.residual_factor),
                materialize_build(build_backend, screen, row_positive.operator.cov_zz),
                materialize_build(build_backend, screen, row_positive.operator.cov_xx),
                materialize_build(build_backend, screen, row_positive.operator.cov_xz),
                materialize_build(build_backend, screen, row_positive.operator.cov_zx),
                materialize_build(build_backend, screen, row_positive.operator.residual_covariance),
                materialize_build(build_backend, Vector{T}(undef, length(row_positive.operator.singular_values)),
                    row_positive.operator.singular_values),
                T(row_positive.operator.condition_ratio),
                row_positive.operator.orientation,
                row_positive.operator.side,
            ),
        ),
        InfiniteBoundaryModel(
            InfiniteBoundaryStencil(
                row_negative.stencil.stencil_coords,
                row_negative.stencil.boundary_coords,
                row_negative.stencil.stencil_positions,
                row_negative.stencil.boundary_positions,
                row_negative.stencil.orientation,
                row_negative.stencil.side,
            ),
            InfiniteBoundaryOperator(
                materialize_build(build_backend, screen, row_negative.operator.predictor),
                materialize_build(build_backend, screen, row_negative.operator.residual_factor),
                materialize_build(build_backend, screen, row_negative.operator.cov_zz),
                materialize_build(build_backend, screen, row_negative.operator.cov_xx),
                materialize_build(build_backend, screen, row_negative.operator.cov_xz),
                materialize_build(build_backend, screen, row_negative.operator.cov_zx),
                materialize_build(build_backend, screen, row_negative.operator.residual_covariance),
                materialize_build(build_backend, Vector{T}(undef, length(row_negative.operator.singular_values)),
                    row_negative.operator.singular_values),
                T(row_negative.operator.condition_ratio),
                row_negative.operator.orientation,
                row_negative.operator.side,
            ),
        ),
        false,
    )
    screen_model = InfinitePhaseScreen(params, state, generator, screen_telescope)
    _warmup_gpu_infinite_screen!(execution_style(state.screen), screen_model, MersenneTwister(0))
    return screen_model
end

function InfiniteMultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    L0::Real=25.0,
    fractional_cn2::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction::AbstractVector,
    altitude::AbstractVector,
    screen_resolution::Int=default_infinite_screen_resolution(tel.params.resolution),
    stencil_size::Int=default_infinite_stencil_size(tel.params.resolution),
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)
    n_layers = length(fractional_cn2)
    n_layers > 0 || throw(InvalidConfiguration("fractional_cn2 cannot be empty"))
    if length(wind_speed) != n_layers || length(wind_direction) != n_layers || length(altitude) != n_layers
        throw(InvalidConfiguration("layer parameter lengths must match fractional_cn2"))
    end
    all(>=(0), fractional_cn2) || throw(InvalidConfiguration("fractional_cn2 must be non-negative"))
    isapprox(sum(fractional_cn2), 1; atol=1e-6, rtol=1e-6) ||
        throw(InvalidConfiguration("fractional_cn2 must sum to 1"))

    params = InfiniteMultiLayerParams(
        T.(fractional_cn2),
        T.(wind_speed),
        T.(wind_direction),
        T.(altitude),
        T[T(wind_speed[i]) * cosd(T(wind_direction[i])) for i in 1:n_layers],
        T[T(wind_speed[i]) * sind(T(wind_direction[i])) for i in 1:n_layers],
        T(r0),
        T(L0),
    )
    layers = [
        InfiniteAtmosphereLayer(
            InfiniteLayerParams(
                T(sqrt(params.cn2_fractions[i])),
                params.wind_velocity_x[i],
                params.wind_velocity_y[i],
                params.altitude[i],
            ),
            InfinitePhaseScreen(tel;
                r0=params.r0,
                L0=params.L0,
                screen_resolution=screen_resolution,
                stencil_size=stencil_size,
                T=T,
                backend=backend,
            ),
            InfiniteLayerState(zero(T), zero(T), 0, 0),
        ) for i in 1:n_layers
    ]
    opd = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    layer_buffer = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    fill!(opd, zero(T))
    fill!(layer_buffer, zero(T))
    state = InfiniteMultiLayerState{T, typeof(opd)}(opd, layer_buffer)
    return InfiniteMultiLayerAtmosphere(params, layers, state)
end

function advance!(atm::InfiniteMultiLayerAtmosphere, tel::Telescope, rng::AbstractRNG)
    fill!(atm.state.opd, zero(eltype(atm.state.opd)))

    for layer in atm.layers
        sample_layer!(atm.state.layer_buffer, layer, tel, rng)
        atm.state.opd .+= atm.state.layer_buffer
    end

    return atm
end

function advance!(atm::InfiniteMultiLayerAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    return advance!(atm, tel, rng)
end

function propagate!(atm::InfiniteMultiLayerAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end

function propagate!(atm::InfiniteMultiLayerAtmosphere, tel::Telescope, src::AbstractSource)
    fill!(atm.state.opd, zero(eltype(atm.state.opd)))
    for layer in atm.layers
        render_layer!(atm.state.layer_buffer, layer, tel, src)
        atm.state.opd .+= atm.state.layer_buffer
    end
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end

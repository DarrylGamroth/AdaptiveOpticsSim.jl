@kernel function fill_telescope_field_kernel!(out, pupil_reflectivity, opd, phase_shift, amp_scale,
    opd_to_cycles, ox::Int, oy::Int, n::Int, n_pad::Int, center_even_grid::Bool)
    i, j = @index(Global, NTuple)
    if i <= n_pad && j <= n_pad
        xi = i - ox
        yj = j - oy
        val = zero(eltype(out))
        if 1 <= xi <= n && 1 <= yj <= n
            @inbounds val = amp_scale * sqrt(pupil_reflectivity[xi, yj]) * cispi(opd_to_cycles * opd[xi, yj])
        end
        if center_even_grid
            val *= cis(phase_shift * (i + j - 2))
        end
        @inbounds out[i, j] = val
    end
end

@kernel function fill_pupil_field_kernel!(out, amplitude, opd, phase_shift,
    amplitude_scale, opd_to_cycles, ox::Int, oy::Int, n::Int, n_pad::Int,
    apply_centering::Bool)
    i, j = @index(Global, NTuple)
    if i <= n_pad && j <= n_pad
        xi = i - ox
        yj = j - oy
        value = zero(eltype(out))
        if 1 <= xi <= n && 1 <= yj <= n
            @inbounds value = amplitude_scale * amplitude[xi, yj] *
                cispi(opd_to_cycles * opd[xi, yj])
        end
        if apply_centering
            value *= cis(phase_shift * (i + j - 2))
        end
        @inbounds out[i, j] = value
    end
end

@kernel function apply_phase_opd_kernel!(field, phase_or_opd, opd_to_cycles, ox::Int, oy::Int, n::Int,
    full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(phase_or_opd, 1) && j <= size(phase_or_opd, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= cispi(opd_to_cycles * phase_or_opd[i, j])
    end
end

@kernel function apply_phase_rad_kernel!(field, phase_or_opd, ox::Int, oy::Int, n::Int, full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(phase_or_opd, 1) && j <= size(phase_or_opd, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= cis(phase_or_opd[i, j])
    end
end

@kernel function apply_amplitude_kernel!(field, amplitude, ox::Int, oy::Int, n::Int, full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(amplitude, 1) && j <= size(amplitude, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= amplitude[i, j]
    end
end

@kernel function intensity_kernel!(out, field, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds out[i, j] = abs2(field[i, j])
    end
end

@kernel function accumulate_abs2_kernel!(out, field, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds out[i, j] += abs2(field[i, j])
    end
end

struct PupilFieldFormationPlan{
    T<:AbstractFloat,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
}
    input_metadata::I
    output_metadata::O
    embedding_offsets::NTuple{2,Int}
    amplitude_scale::T
    opd_to_cycles::T
    centering_phase_shift::T
    apply_centering_phase::Bool
end

@inline function field_embedding_offsets(resolution::Int, padded_resolution::Int)
    return div(padded_resolution - resolution, 2), div(padded_resolution - resolution, 2)
end

@inline function field_active_axes(plan::PupilFieldFormationPlan)
    ox, oy = plan.embedding_offsets
    n, m = plan.input_metadata.dimensions
    return (ox + 1:ox + n, oy + 1:oy + m)
end

@inline source_field_normalization(::PhysicalPhotonIrradianceSource) =
    PhotonRateNormalization()
@inline source_field_normalization(::NormalizedTestSource) =
    DimensionlessNormalization()
@inline source_field_measure(::AbstractSourceRadiometry) =
    CellIntegratedMeasure()

function ElectricField(wavefront::PupilFunction, src::AbstractSource;
    zero_padding::Int=1,
    T::Type{<:AbstractFloat}=eltype(wavefront.opd),
    normalization::AbstractOpticalNormalization=
        source_field_normalization(source_radiometry(src)),
    spatial_measure::AbstractSpatialMeasure=
        source_field_measure(source_radiometry(src)),
    coherence::AbstractCombinationPolicy=CoherentFieldCombination())
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n, m = wavefront.metadata.dimensions
    n == m || throw(DimensionMismatchError(
        "PupilFunction must be square to prepare an ElectricField"))
    n_pad = n * zero_padding
    values = similar(wavefront.opd, Complex{T}, n_pad, n_pad)
    fill!(values, zero(eltype(values)))
    sampling = (T(wavefront.metadata.sampling[1]),
        T(wavefront.metadata.sampling[2]))
    metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=MetricCoordinates(),
        sampling=sampling,
        orientation=wavefront.metadata.orientation,
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=normalization,
        spatial_measure=spatial_measure,
        coherence=coherence)
    return ElectricField(metadata, values)
end

@inline electric_field_wavelength(field::ElectricField) =
    electric_field_wavelength(field.metadata.spectral)
@inline electric_field_wavelength(channel::MonochromaticChannel) =
    channel.wavelength_m
function electric_field_wavelength(::AbstractSpectralCoordinate)
    throw(InvalidConfiguration(
        "ElectricField must declare a monochromatic wavelength"))
end

function prepare_pupil_field(tel::Telescope, wavefront::PupilFunction,
    src::AbstractSource, field::ElectricField;
    center_even_grid::Bool=true,
    amplitude_scale::Union{Real,Nothing}=nothing)
    validate_plane_storage(wavefront.metadata, wavefront.amplitude;
        label="pupil-wavefront amplitude")
    validate_plane_storage(wavefront.metadata, wavefront.opd;
        label="pupil-wavefront OPD")
    validate_plane_storage(field.metadata, field.values;
        label="electric field")
    require_centered_plane_geometry(wavefront.metadata;
        label="PupilFunction")
    require_centered_plane_geometry(field.metadata;
        label="ElectricField")
    require_metric_coordinates(wavefront.metadata;
        label="pupil-field input PupilFunction")
    require_metric_coordinates(field.metadata;
        label="pupil-field output ElectricField")
    typeof(wavefront.metadata.kind) === PupilPlane ||
        throw(InvalidConfiguration(
            "pupil-field formation requires a pupil-plane wavefront"))
    typeof(field.metadata.kind) === PupilPlane ||
        throw(InvalidConfiguration(
            "pupil-field formation requires a pupil-plane ElectricField"))
    wavefront.metadata.dimensions == size(pupil_mask(tel)) ||
        throw(DimensionMismatchError(
            "PupilFunction dimensions must match telescope aperture"))
    plane_device(pupil_mask(tel)) == wavefront.metadata.device ||
        throw(InvalidConfiguration(
            "telescope aperture and PupilFunction occupy different physical devices"))
    typeof(backend(tel)) === typeof(backend(wavefront)) ===
        typeof(backend(field)) || throw(InvalidConfiguration(
            "telescope, PupilFunction, and ElectricField backends must match"))
    wavefront.metadata.device == field.metadata.device ||
        throw(InvalidConfiguration(
            "PupilFunction and ElectricField must occupy the same physical device"))
    wavefront.metadata.sampling == field.metadata.sampling ||
        throw(InvalidConfiguration(
            "PupilFunction and ElectricField sampling must match"))
    wavefront.metadata.orientation == field.metadata.orientation ||
        throw(InvalidConfiguration(
            "PupilFunction and ElectricField axis orientation must match"))
    eltype(wavefront.opd) === real(eltype(field.values)) ||
        throw(InvalidConfiguration(
            "PupilFunction and ElectricField real numeric types must match"))
    field.metadata.spectral == MonochromaticChannel(
        eltype(wavefront.opd)(wavelength(src))) ||
        throw(InvalidConfiguration(
            "source wavelength must match ElectricField wavelength"))

    n, m = wavefront.metadata.dimensions
    n_pad, m_pad = field.metadata.dimensions
    n == m && n_pad == m_pad || throw(DimensionMismatchError(
        "pupil-field formation requires square input and output grids"))
    n_pad % n == 0 || throw(DimensionMismatchError(
        "ElectricField dimensions must be an integer multiple of the pupil grid"))
    T = eltype(wavefront.opd)
    ox, oy = field_embedding_offsets(n, n_pad)
    apply_centering = center_even_grid && iseven(n_pad)
    phase_shift = apply_centering ?
        -T(pi) * (T(n_pad) + one(T)) / T(n_pad) : zero(T)
    pixel_area = wavefront.metadata.sampling[1] *
        wavefront.metadata.sampling[2]
    resolved_amplitude_scale = _resolve_pupil_amplitude_scale(
        source_radiometry(src), wavefront, field.metadata, src,
        pixel_area, amplitude_scale, T)
    resolved_amplitude_scale >= zero(T) || throw(InvalidConfiguration(
        "pupil-field amplitude scale must be non-negative"))
    return PupilFieldFormationPlan{
        T,typeof(wavefront.metadata),typeof(field.metadata),
    }(
        wavefront.metadata,
        field.metadata,
        (ox, oy),
        resolved_amplitude_scale,
        T(2) / T(wavelength(src)),
        phase_shift,
        apply_centering,
    )
end

function _resolve_pupil_amplitude_scale(
    radiometry::AbstractSourceRadiometry, wavefront::PupilFunction,
    metadata::OpticalPlaneMetadata, src::AbstractSource, pixel_area,
    ::Nothing, ::Type{T}) where {T<:AbstractFloat}
    _require_source_field_contract(radiometry, metadata)
    return _default_pupil_amplitude_scale(radiometry, wavefront, src,
        pixel_area, T)
end

function _resolve_pupil_amplitude_scale(::AbstractSourceRadiometry,
    ::PupilFunction, metadata::OpticalPlaneMetadata, ::AbstractSource,
    pixel_area, amplitude_scale::Real, ::Type{T}) where {T<:AbstractFloat}
    _require_dimensionless_field_normalization(metadata.normalization)
    _require_dimensionless_field_measure(metadata.spatial_measure)
    _require_coherent_field(metadata.coherence)
    return T(amplitude_scale)
end

@inline _require_dimensionless_field_normalization(
    ::DimensionlessNormalization) = nothing
function _require_dimensionless_field_normalization(
    ::AbstractOpticalNormalization)
    throw(InvalidConfiguration(
        "an explicit pupil-field amplitude scale requires dimensionless field metadata"))
end

@inline _require_dimensionless_field_measure(::PointSampledMeasure) = nothing
@inline _require_dimensionless_field_measure(::CellIntegratedMeasure) = nothing
function _require_dimensionless_field_measure(::AbstractSpatialMeasure)
    throw(InvalidConfiguration(
        "dimensionless source fields must be point-sampled or cell-integrated"))
end


@inline function _require_source_field_contract(
    ::PhysicalPhotonIrradianceSource, metadata::OpticalPlaneMetadata)
    _require_physical_field_normalization(metadata.normalization)
    _require_cell_integrated_field(metadata.spatial_measure)
    _require_coherent_field(metadata.coherence)
    return nothing
end

@inline function _require_source_field_contract(::NormalizedTestSource,
    metadata::OpticalPlaneMetadata)
    _require_normalized_field_normalization(metadata.normalization)
    _require_cell_integrated_field(metadata.spatial_measure)
    _require_coherent_field(metadata.coherence)
    return nothing
end

@inline _require_physical_field_normalization(::PhotonRateNormalization) =
    nothing
function _require_physical_field_normalization(::AbstractOpticalNormalization)
    throw(InvalidConfiguration(
        "physical sources require photon-rate field normalization"))
end

@inline _require_normalized_field_normalization(
    ::DimensionlessNormalization) = nothing
function _require_normalized_field_normalization(
    ::AbstractOpticalNormalization)
    throw(InvalidConfiguration(
        "normalized sources require dimensionless field normalization"))
end

@inline _require_cell_integrated_field(::CellIntegratedMeasure) = nothing
function _require_cell_integrated_field(::AbstractSpatialMeasure)
    throw(InvalidConfiguration(
        "source fields must use a cell-integrated spatial measure"))
end

@inline _require_coherent_field(::CoherentFieldCombination) = nothing
function _require_coherent_field(::AbstractCombinationPolicy)
    throw(InvalidConfiguration(
        "source fields must declare coherent field combination"))
end


@inline function _default_pupil_amplitude_scale(
    ::PhysicalPhotonIrradianceSource, ::PupilFunction, src::AbstractSource,
    pixel_area, ::Type{T}) where {T<:AbstractFloat}
    return sqrt(T(photon_irradiance(src) * pixel_area))
end


function _default_pupil_amplitude_scale(::NormalizedTestSource,
    wavefront::PupilFunction, src::AbstractSource, pixel_area,
    ::Type{T}) where {T<:AbstractFloat}
    transmitted = T(sum(abs2, wavefront.amplitude))
    transmitted > zero(T) || throw(InvalidConfiguration(
        "normalized pupil field requires non-zero transmitting amplitude"))
    return sqrt(T(source_radiometric_value(src)) / transmitted)
end

function fill_telescope_field!(out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource;
    zero_padding::Int=1,
    center_even_grid::Bool=true) where {T<:AbstractFloat}
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n = tel.params.resolution
    n_pad = n * zero_padding
    size(out) == (n_pad, n_pad) ||
        throw(DimensionMismatchError("electric field size must match telescope resolution * zero_padding"))

    _fill_telescope_field!(execution_style(out), out, tel, src, zero_padding, center_even_grid)
    return out
end

function fill_telescope_field_async!(out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource;
    zero_padding::Int=1,
    center_even_grid::Bool=true) where {T<:AbstractFloat}
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n = tel.params.resolution
    n_pad = n * zero_padding
    size(out) == (n_pad, n_pad) ||
        throw(DimensionMismatchError("electric field size must match telescope resolution * zero_padding"))
    _fill_telescope_field_async!(execution_style(out), out, tel, src, zero_padding, center_even_grid)
    return out
end

function _fill_telescope_field!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    fill!(out, zero(eltype(out)))
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    reflectivity = pupil_reflectivity(tel)
    @views @. out[ox+1:ox+n, oy+1:oy+n] = amp_scale *
        sqrt(reflectivity) * cispi(opd_to_cycles * tel.state.opd)
    if center_even_grid && iseven(n_pad)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(ScalarCPUStyle(), out, phase_shift)
    end
    return out
end

function _fill_telescope_field!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    phase_shift = center_even_grid && iseven(n_pad) ? -T(pi) * (T(n_pad) + one(T)) / T(n_pad) : zero(T)
    launch_kernel!(style, fill_telescope_field_kernel!, out, pupil_reflectivity(tel), tel.state.opd, phase_shift,
        amp_scale, opd_to_cycles, ox, oy, n, n_pad, center_even_grid && iseven(n_pad); ndrange=size(out))
    return out
end

function _fill_telescope_field_async!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    return _fill_telescope_field!(ScalarCPUStyle(), out, tel, src, zero_padding, center_even_grid)
end

function _fill_telescope_field_async!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    phase_shift = center_even_grid && iseven(n_pad) ? -T(pi) * (T(n_pad) + one(T)) / T(n_pad) : zero(T)
    launch_kernel_async!(style, fill_telescope_field_kernel!, out, pupil_reflectivity(tel), tel.state.opd, phase_shift,
        amp_scale, opd_to_cycles, ox, oy, n, n_pad, center_even_grid && iseven(n_pad); ndrange=size(out))
    return out
end

function fill_electric_field!(field::ElectricField,
    wavefront::PupilFunction, plan::PupilFieldFormationPlan)
    field.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "ElectricField metadata does not match its prepared formation plan"))
    wavefront.metadata == plan.input_metadata || throw(InvalidConfiguration(
        "PupilFunction metadata does not match its prepared formation plan"))
    _fill_electric_field!(execution_style(field.values), field, wavefront,
        plan)
    return field
end

function fill_electric_field_async!(field::ElectricField,
    wavefront::PupilFunction, plan::PupilFieldFormationPlan)
    field.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "ElectricField metadata does not match its prepared formation plan"))
    wavefront.metadata == plan.input_metadata || throw(InvalidConfiguration(
        "PupilFunction metadata does not match its prepared formation plan"))
    _fill_electric_field_async!(execution_style(field.values), field,
        wavefront, plan)
    return field
end

function _fill_electric_field!(::ScalarCPUStyle, field::ElectricField,
    wavefront::PupilFunction, plan::PupilFieldFormationPlan)
    out = field.values
    fill!(out, zero(eltype(out)))
    axes = field_active_axes(plan)
    @views @. out[axes...] = plan.amplitude_scale * wavefront.amplitude *
        cispi(plan.opd_to_cycles * wavefront.opd)
    if plan.apply_centering_phase
        apply_centering_phase!(ScalarCPUStyle(), out,
            plan.centering_phase_shift)
    end
    return field
end


function _fill_electric_field!(style::AcceleratorStyle,
    field::ElectricField, wavefront::PupilFunction,
    plan::PupilFieldFormationPlan)
    n = plan.input_metadata.dimensions[1]
    n_pad = plan.output_metadata.dimensions[1]
    ox, oy = plan.embedding_offsets
    launch_kernel!(style, fill_pupil_field_kernel!, field.values,
        wavefront.amplitude, wavefront.opd, plan.centering_phase_shift,
        plan.amplitude_scale, plan.opd_to_cycles, ox, oy, n, n_pad,
        plan.apply_centering_phase; ndrange=size(field.values))
    return field
end

function _fill_electric_field_async!(::ScalarCPUStyle,
    field::ElectricField, wavefront::PupilFunction,
    plan::PupilFieldFormationPlan)
    return _fill_electric_field!(ScalarCPUStyle(), field, wavefront, plan)
end

function _fill_electric_field_async!(style::AcceleratorStyle,
    field::ElectricField, wavefront::PupilFunction,
    plan::PupilFieldFormationPlan)
    n = plan.input_metadata.dimensions[1]
    n_pad = plan.output_metadata.dimensions[1]
    ox, oy = plan.embedding_offsets
    launch_kernel_async!(style, fill_pupil_field_kernel!, field.values,
        wavefront.amplitude, wavefront.opd, plan.centering_phase_shift,
        plan.amplitude_scale, plan.opd_to_cycles, ox, oy, n, n_pad,
        plan.apply_centering_phase; ndrange=size(field.values))
    return field
end

function field_target_view(field::ElectricField, input::AbstractMatrix,
    plan::Union{PupilFieldFormationPlan,Nothing}=nothing)
    if size(input) == size(field.values)
        return field.values
    elseif !isnothing(plan) && size(input) == plan.input_metadata.dimensions
        active_axes = field_active_axes(plan)
        return @view field.values[active_axes...]
    end
    throw(DimensionMismatchError(
        "field map size must match the full field or prepared active pupil"))
end

function apply_phase!(field::ElectricField, phase_or_opd::AbstractMatrix; units::Symbol=:opd)
    _apply_phase!(execution_style(field.values), field, phase_or_opd, units,
        nothing)
    return field
end

function apply_phase!(field::ElectricField, phase_or_opd::AbstractMatrix,
    plan::PupilFieldFormationPlan; units::Symbol=:opd)
    _apply_phase!(execution_style(field.values), field, phase_or_opd, units,
        plan)
    return field
end

function apply_phase_async!(field::ElectricField,
    phase_or_opd::AbstractMatrix,
    plan::Union{PupilFieldFormationPlan,Nothing}=nothing; units::Symbol=:opd)
    _apply_phase_async!(execution_style(field.values), field, phase_or_opd,
        units, plan)
    return field
end

function _phase_target_layout(field::ElectricField,
    input::AbstractMatrix,
    plan::Union{PupilFieldFormationPlan,Nothing})
    if size(input) == size(field.values)
        return true, 0, 0
    elseif !isnothing(plan) && size(input) == plan.input_metadata.dimensions
        ox, oy = plan.embedding_offsets
        return false, ox, oy
    end
    throw(DimensionMismatchError(
        "field map size must match the full field or prepared active pupil"))
end

function _apply_phase!(::ScalarCPUStyle, field::ElectricField,
    phase_or_opd::AbstractMatrix, units::Symbol,
    plan::Union{PupilFieldFormationPlan,Nothing})
    target = field_target_view(field, phase_or_opd, plan)
    if units === :opd
        T = real(eltype(field.values))
        opd_to_cycles = T(2) / T(electric_field_wavelength(field))
        @. target *= cispi(opd_to_cycles * phase_or_opd)
        return field
    elseif units === :phase
        @. target *= cis(phase_or_opd)
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function _apply_phase!(style::AcceleratorStyle, field::ElectricField,
    phase_or_opd::AbstractMatrix, units::Symbol,
    plan::Union{PupilFieldFormationPlan,Nothing})
    full_field, ox, oy = _phase_target_layout(field, phase_or_opd, plan)
    if units === :opd
        T = real(eltype(field.values))
        opd_to_cycles = T(2) / T(electric_field_wavelength(field))
        launch_kernel!(style, apply_phase_opd_kernel!, field.values,
            phase_or_opd, opd_to_cycles, ox, oy,
            isnothing(plan) ? size(field.values, 1) :
                plan.input_metadata.dimensions[1],
            full_field; ndrange=size(phase_or_opd))
        return field
    elseif units === :phase
        launch_kernel!(style, apply_phase_rad_kernel!, field.values,
            phase_or_opd, ox, oy,
            isnothing(plan) ? size(field.values, 1) :
                plan.input_metadata.dimensions[1],
            full_field; ndrange=size(phase_or_opd))
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function _apply_phase_async!(::ScalarCPUStyle, field::ElectricField,
    phase_or_opd::AbstractMatrix, units::Symbol,
    plan::Union{PupilFieldFormationPlan,Nothing})
    return _apply_phase!(ScalarCPUStyle(), field, phase_or_opd, units, plan)
end

function _apply_phase_async!(style::AcceleratorStyle,
    field::ElectricField, phase_or_opd::AbstractMatrix, units::Symbol,
    plan::Union{PupilFieldFormationPlan,Nothing})
    full_field, ox, oy = _phase_target_layout(field, phase_or_opd, plan)
    if units === :opd
        T = real(eltype(field.values))
        opd_to_cycles = T(2) / T(electric_field_wavelength(field))
        launch_kernel_async!(style, apply_phase_opd_kernel!, field.values,
            phase_or_opd, opd_to_cycles, ox, oy,
            isnothing(plan) ? size(field.values, 1) :
                plan.input_metadata.dimensions[1],
            full_field; ndrange=size(phase_or_opd))
        return field
    elseif units === :phase
        launch_kernel_async!(style, apply_phase_rad_kernel!, field.values,
            phase_or_opd, ox, oy,
            isnothing(plan) ? size(field.values, 1) :
                plan.input_metadata.dimensions[1],
            full_field; ndrange=size(phase_or_opd))
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function apply_amplitude!(field::ElectricField, amplitude::AbstractMatrix)
    _apply_amplitude!(execution_style(field.values), field, amplitude, nothing)
    return field
end

function apply_amplitude!(field::ElectricField, amplitude::AbstractMatrix,
    plan::PupilFieldFormationPlan)
    _apply_amplitude!(execution_style(field.values), field, amplitude, plan)
    return field
end

function _apply_amplitude!(::ScalarCPUStyle, field::ElectricField,
    amplitude::AbstractMatrix,
    plan::Union{PupilFieldFormationPlan,Nothing})
    target = field_target_view(field, amplitude, plan)
    @. target *= amplitude
    return field
end

function _apply_amplitude!(style::AcceleratorStyle, field::ElectricField,
    amplitude::AbstractMatrix,
    plan::Union{PupilFieldFormationPlan,Nothing})
    full_field, ox, oy = _phase_target_layout(field, amplitude, plan)
    launch_kernel!(style, apply_amplitude_kernel!, field.values, amplitude,
        ox, oy,
        isnothing(plan) ? size(field.values, 1) :
            plan.input_metadata.dimensions[1],
        full_field; ndrange=size(amplitude))
    return field
end

function intensity!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    require_same_backend(out, field)
    _intensity!(execution_style(out), out, field.values)
    return out
end

function intensity_async!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    require_same_backend(out, field)
    _intensity_async!(execution_style(out), out, field.values)
    return out
end

function _intensity!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out = abs2(field_values)
    return out
end

function _intensity!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel!(style, intensity_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function _intensity_async!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out = abs2(field_values)
    return out
end

function _intensity_async!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel_async!(style, intensity_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function _accumulate_field_intensity!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    require_same_backend(out, field)
    _accumulate_intensity!(execution_style(out), out, field.values)
    return out
end

function _accumulate_field_intensity_async!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.values) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    require_same_backend(out, field)
    _accumulate_intensity_async!(execution_style(out), out, field.values)
    return out
end

function _accumulate_intensity!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out += abs2(field_values)
    return out
end

function _accumulate_intensity!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel!(style, accumulate_abs2_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function _accumulate_intensity_async!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out += abs2(field_values)
    return out
end

function _accumulate_intensity_async!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel_async!(style, accumulate_abs2_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

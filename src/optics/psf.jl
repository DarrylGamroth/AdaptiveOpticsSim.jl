@kernel function apply_centering_phase_kernel!(field, phase_shift, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds field[i, j] *= cis(phase_shift * (i + j - 2))
    end
end

function apply_centering_phase!(::ScalarCPUStyle, field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    n, m = size(field)
    @inbounds for j in 1:m, i in 1:n
        field[i, j] *= cis(phase_shift * (i + j - 2))
    end
    return field
end

function apply_centering_phase!(style::AcceleratorStyle, field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    launch_kernel!(style, apply_centering_phase_kernel!, field, phase_shift, size(field, 1); ndrange=size(field))
    return field
end

function ensure_psf_state!(tel::Telescope, n::Int)
    if size(tel.state.psf) != (n, n)
        tel.state.psf = similar(tel.state.psf, n, n)
    end
    return tel
end

function ensure_psf_workspace!(tel::Telescope, n::Int)
    ensure_psf_buffers!(tel.state.psf_workspace, n)
    return tel.state.psf_workspace
end

function centered_psf_from_field!(out::AbstractMatrix{T},
    field::ElectricField,
    propagation::FraunhoferPropagation) where {T<:AbstractFloat}
    return fraunhofer_intensity_from_field!(out, field, propagation)
end

struct DirectPSFPlan{
    T<:AbstractFloat,
    F<:PupilFieldFormationPlan,
    O<:OpticalPlaneMetadata,
}
    field_formation::F
    output_metadata::O
    shift_pixels::NTuple{2,T}
end

struct DirectPSFWorkspace{
    P<:FraunhoferPropagation,
    R<:AbstractMatrix,
}
    propagation::P
    unshifted_intensity::R
end

function prepare_direct_psf(tel::Telescope, wavefront::PupilFunction,
    src::Source, field::ElectricField, output::IntensityMap)
    propagation = FraunhoferPropagation(field)
    return _prepare_direct_psf(tel, wavefront, src, field, output,
        propagation)
end

function prepare_direct_psf(tel::Telescope, wavefront::PupilFunction,
    src::Source, field::ElectricField)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    prepared = _prepare_direct_psf(tel, wavefront, src, field, output,
        propagation)
    return (; output, prepared.plan, prepared.workspace)
end

function _prepare_direct_psf(tel::Telescope, wavefront::PupilFunction,
    src::Source, field::ElectricField, output::IntensityMap,
    propagation::FraunhoferPropagation)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    require_same_plane_grid(output.metadata, propagation.output_metadata;
        label="direct-PSF output", require_numeric_type=false)
    require_compatible_radiometry(output.metadata,
        propagation.output_metadata; label="direct-PSF output")
    _require_incoherent_policy(output.metadata.coherence,
        "direct-PSF output")
    scale_arcsec = propagation.params.output_sampling_rad *
        (180 / pi) * 3600
    coordinates = src.params.coordinates_xy_arcsec
    T = eltype(output.values)
    shift_pixels = (
        T(coordinates[1] / scale_arcsec),
        T(coordinates[2] / scale_arcsec),
    )
    unshifted = similar(output.values)
    plan = DirectPSFPlan{
        T,typeof(formation),typeof(output.metadata),
    }(formation, output.metadata, shift_pixels)
    workspace = DirectPSFWorkspace(propagation, unshifted)
    return (; plan, workspace)
end

function compute_psf!(output::IntensityMap, field::ElectricField,
    wavefront::PupilFunction, plan::DirectPSFPlan,
    workspace::DirectPSFWorkspace)
    output.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "IntensityMap metadata does not match its prepared direct-PSF plan"))
    fill_electric_field!(field, wavefront, plan.field_formation)
    dx, dy = plan.shift_pixels
    if iszero(dx) && iszero(dy)
        fraunhofer_intensity_from_field!(output.values, field,
            workspace.propagation)
    else
        fraunhofer_intensity_from_field!(workspace.unshifted_intensity,
            field, workspace.propagation)
        shift_psf!(output.values, workspace.unshifted_intensity, dx, dy)
    end
    return output
end

function compute_psf_centered!(tel::Telescope, src::Source, ws::Workspace, zero_padding::Int=1)
    _require_physical_photon_irradiance(src, "legacy telescope PSF")
    n = tel.params.resolution
    if zero_padding < 1
        throw(InvalidConfiguration("zero_padding must be >= 1"))
    end
    n_pad = n * zero_padding
    T = eltype(tel.state.opd)

    ensure_psf_buffers!(ws, n_pad)
    fill_telescope_field!(ws.pupil_field, tel, src; zero_padding=zero_padding)

    copyto!(ws.fft_buffer, ws.pupil_field)
    execute_fft_plan!(ws.fft_buffer, ws.fft_plan)

    fft_scale = inv(T(n_pad))
    @. ws.psf_buffer = abs2(ws.fft_buffer) * (fft_scale * fft_scale)

    ensure_psf_state!(tel, n_pad)
    copyto!(tel.state.psf, ws.psf_buffer)
    return tel.state.psf
end

function compute_psf!(tel::Telescope, src::Source, ws::Workspace, zero_padding::Int=1)
    psf = compute_psf_centered!(tel, src, ws, zero_padding)
    coords_xy_arcsec = src.params.coordinates_xy_arcsec
    if iszero(coords_xy_arcsec[1]) && iszero(coords_xy_arcsec[2])
        return psf
    end
    scale = psf_pixel_scale_arcsec(tel, src, zero_padding)
    shift_psf!(tel.state.psf, ws.psf_buffer, coords_xy_arcsec[1] / scale, coords_xy_arcsec[2] / scale)
    return tel.state.psf
end

function compute_psf!(tel::Telescope, src::Source, zero_padding::Int)
    _require_physical_photon_irradiance(src, "legacy telescope PSF")
    if zero_padding < 1
        throw(InvalidConfiguration("zero_padding must be >= 1"))
    end
    n_pad = tel.params.resolution * zero_padding
    ws = ensure_psf_workspace!(tel, n_pad)
    return compute_psf!(tel, src, ws, zero_padding)
end

function compute_psf!(tel::Telescope, src::Source; zero_padding::Int=1, ws::Union{Workspace,Nothing}=nothing)
    if ws === nothing
        return compute_psf!(tel, src, zero_padding)
    end
    return compute_psf!(tel, src, ws, zero_padding)
end

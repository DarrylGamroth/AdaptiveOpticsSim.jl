import ForwardDiff

const MISREGISTRATION_AD_FIELDS = (
    :shift_x,
    :shift_y,
    :rotation_deg,
    :radial_scaling,
    :tangential_scaling,
)

function _compute_meta_sensitivity_matrix_ad(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS,
    basis::AbstractMatrix; source=nothing,
    misregistration_zero::Misregistration=Misregistration(T=eltype(tel.state.opd)),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(tel.state.opd)),
    direction_epsilon::Real=sqrt(eps(eltype(tel.state.opd))),
    n_mis_reg::Int=3, field_order=collect(MISREG_FIELDS),
    amplitude::Real=1e-9, cache_path::Union{Nothing,String}=nothing,
    save_sensitivity::Bool=true, recompute_sensitivity::Bool=false,
    wfs_mis_registered::Bool=false)

    if cache_path !== nothing && !recompute_sensitivity && isfile(cache_path)
        return deserialize(cache_path)
    end
    wfs_mis_registered &&
        throw(UnsupportedAlgorithm(
            "AD sensitivity supports DM misregistration only; use sensitivity=:finite_difference for WFS misregistration"))
    _require_cpu_ad_probe(tel, dm, wfs)

    T = eltype(tel.state.opd)
    fields = Tuple(collect(field_order)[1:min(n_mis_reg, length(field_order))])
    _validate_misregistration_ad_field_order(fields)

    dm_model = influence_model(dm)
    dm_topology = topology(dm)
    supports_dm_misregistration_identification(dm_model, dm_topology) ||
        throw(UnsupportedAlgorithm(
            "misregistration AD sensitivity is only supported for grid-backed Gaussian DeformableMirror models"))
    size(basis, 1) == length(dm.state.coefs) ||
        throw(DimensionMismatchError("basis row count must match DM coefficients"))
    direction_epsilon > 0 || throw(InvalidConfiguration("direction_epsilon must be positive"))

    dm0 = DeformableMirror(tel; topology=dm_topology, influence_model=dm_model,
        misregistration=misregistration_zero, T=T)
    calib0 = _interaction_matrix_for_sensitivity(dm0, wfs, tel, basis, source, amplitude)
    calib0_control_matrix = ControlMatrix(calib0.matrix)

    d_modes = _gaussian_dm_mode_parameter_jacobians(tel, dm0, fields)
    base_modes = sampled_influence_matrix(dm0)
    meta = _wfs_directional_meta_sensitivity(tel, wfs, basis, source, calib0.matrix,
        base_modes, d_modes, T(amplitude), T(direction_epsilon))

    meta_control_matrix = ControlMatrix(meta)
    out = MetaSensitivity(meta_control_matrix, calib0_control_matrix, epsilon, collect(fields))
    if cache_path !== nothing && save_sensitivity
        serialize(cache_path, out)
    end
    return out
end

function compute_meta_sensitivity_matrix_ad_probe(args...; kwargs...)
    return compute_meta_sensitivity_matrix(args...; sensitivity=:ad, kwargs...)
end

function _require_cpu_ad_probe(tel, dm, wfs=nothing)
    _is_cpu_array(tel.state.opd) && _is_cpu_array(tel.state.pupil) && _is_cpu_array(dm.state.modes) ||
        throw(UnsupportedAlgorithm(
            "ForwardDiff AD sensitivity is CPU-only; use sensitivity=:finite_difference for accelerator-backed arrays"))
    if wfs !== nothing
        _is_cpu_array(slopes(wfs)) ||
            throw(UnsupportedAlgorithm(
                "ForwardDiff AD sensitivity is CPU-only; use sensitivity=:finite_difference for accelerator-backed WFS arrays"))
    end
    return nothing
end

@inline _is_cpu_array(x) = _is_cpu_execution_style(execution_style(x))
@inline _is_cpu_execution_style(::ScalarCPUStyle) = true
@inline _is_cpu_execution_style(::ExecutionStyle) = false

_normalize_ad_field_order(field_order::Symbol) = (field_order,)
_normalize_ad_field_order(field_order::Tuple{Vararg{Symbol}}) = field_order
_normalize_ad_field_order(field_order) = Tuple(Symbol.(field_order))

function _has_duplicate_ad_fields(fields::Tuple)
    @inbounds for i in eachindex(fields)
        field = fields[i]
        for j in 1:(i - 1)
            fields[j] === field && return true
        end
    end
    return false
end

function _validate_misregistration_ad_field_order(fields::Tuple)
    isempty(fields) && throw(InvalidConfiguration("field_order must contain at least one field"))
    !_has_duplicate_ad_fields(fields) ||
        throw(InvalidConfiguration("field_order must not contain duplicates"))
    for field in fields
        field in MISREGISTRATION_AD_FIELDS ||
            throw(InvalidConfiguration("unsupported misregistration AD field $(field)"))
    end
    return fields
end

function _interaction_matrix_for_sensitivity(dm, wfs, tel, basis, source, amplitude)
    if source === nothing
        return interaction_matrix(dm, wfs, tel, basis; amplitude=amplitude)
    end
    return interaction_matrix(dm, wfs, tel, basis, source; amplitude=amplitude)
end

function _measure_for_sensitivity!(wfs, tel, source)
    if source === nothing
        return _measure_for_calibration!(wfs, tel, nothing)
    end
    return _measure_for_calibration!(wfs, tel, source)
end

function _gaussian_dm_mode_parameter_jacobians(tel, dm, fields::Tuple)
    T = eltype(dm.state.modes)
    n_elements, n_commands = size(dm.state.modes)
    out = [Matrix{T}(undef, n_elements, n_commands) for _ in fields]
    base = Vector{T}(undef, length(fields))
    values = Vector{T}(undef, n_elements)
    jacobian = Matrix{T}(undef, n_elements, length(fields))
    cfg = _gaussian_influence_jacobian_config(values, tel, dm, 1, fields, base)
    @inbounds for actuator_index in 1:n_commands
        _base_parameter_vector!(base, tel, dm, actuator_index, fields)
        _gaussian_influence_jacobian!(values, jacobian, tel, dm, actuator_index, fields, base, cfg)
        for field_index in eachindex(fields)
            out[field_index][:, actuator_index] .= @view jacobian[:, field_index]
        end
    end
    return out
end

function _wfs_directional_meta_sensitivity(tel, wfs, basis, source, calib0::AbstractMatrix{T},
    modes::AbstractMatrix{T}, d_modes::AbstractVector, amplitude::T,
    direction_epsilon::T) where {T<:AbstractFloat}
    n_slopes, n_modes = size(calib0)
    n_elements = size(modes, 1)
    meta = zeros(T, length(calib0), length(d_modes))
    base_opd = similar(tel.state.opd)
    d_opd = similar(tel.state.opd)
    base_vec = reshape(base_opd, :)
    d_vec = reshape(d_opd, :)
    slopes_p = Vector{T}(undef, n_slopes)
    slopes_m = Vector{T}(undef, n_slopes)
    saved_opd = copy(tel.state.opd)
    try
        @inbounds for mode_index in 1:n_modes
            command = @view basis[:, mode_index]
            mul!(base_vec, modes, command)
            base_vec .*= amplitude
            for field_index in eachindex(d_modes)
                mul!(d_vec, d_modes[field_index], command)
                d_vec .*= amplitude

                @. tel.state.opd = base_opd + direction_epsilon * d_opd
                _measure_for_sensitivity!(wfs, tel, source)
                copyto!(slopes_p, slopes(wfs))

                @. tel.state.opd = base_opd - direction_epsilon * d_opd
                _measure_for_sensitivity!(wfs, tel, source)
                copyto!(slopes_m, slopes(wfs))

                offset = (mode_index - 1) * n_slopes
                @views meta[offset+1:offset+n_slopes, field_index] .=
                    (slopes_p .- slopes_m) ./ (2 * direction_epsilon)
            end
        end
    finally
        copyto!(tel.state.opd, saved_opd)
    end
    n_elements == length(base_vec) || throw(DimensionMismatchError("DM modes do not match telescope OPD"))
    return meta
end

function _base_parameter_vector!(params::AbstractVector, tel, dm, actuator_index::Int, fields::Tuple)
    T = eltype(params)
    @inbounds for (idx, field) in pairs(fields)
        params[idx] = T(_base_parameter_value(tel, dm, actuator_index, field))
    end
    return params
end

function _base_parameter_value(tel, dm, actuator_index::Int, field::Symbol)
    mis = dm.params.misregistration
    coords = actuator_coordinates(topology(dm))
    if field === :influence_width
        return influence_width(dm)
    elseif field === :mechanical_coupling
        return mechanical_coupling(dm)
    elseif field === :actuator_x
        return coords[1, actuator_index]
    elseif field === :actuator_y
        return coords[2, actuator_index]
    elseif field === :shift_x
        return mis.shift_x
    elseif field === :shift_y
        return mis.shift_y
    elseif field === :rotation_deg
        return rotation_deg(mis)
    elseif field === :radial_scaling
        return mis.radial_scaling
    elseif field === :tangential_scaling
        return mis.tangential_scaling
    end
    throw(InvalidConfiguration("unsupported Gaussian DM AD field $(field)"))
end

@inline function _field_value(fields::Tuple, p::AbstractVector, field::Symbol, default)
    idx = findfirst(==(field), fields)
    return idx === nothing ? convert(eltype(p), default) : p[idx]
end

struct GaussianInfluenceADFunction{TEL,DM,F}
    tel::TEL
    dm::DM
    actuator_index::Int
    fields::F
end

@inline function (f::GaussianInfluenceADFunction)(out::AbstractVector, p::AbstractVector)
    return _gaussian_influence_vector!(out, f.tel, f.dm, f.actuator_index, f.fields, p)
end

function _gaussian_influence_jacobian_config(values, tel, dm, actuator_index::Int, fields::Tuple, base::AbstractVector)
    f = GaussianInfluenceADFunction(tel, dm, actuator_index, fields)
    return ForwardDiff.JacobianConfig(f, values, base)
end

function _gaussian_influence_jacobian!(values::AbstractVector, jacobian::AbstractMatrix, tel, dm,
    actuator_index::Int, fields::Tuple, base::AbstractVector,
    cfg=_gaussian_influence_jacobian_config(values, tel, dm, actuator_index, fields, base))
    f = GaussianInfluenceADFunction(tel, dm, actuator_index, fields)
    ForwardDiff.jacobian!(jacobian, f, values, base, cfg, Val(false))
    return values, jacobian
end

function _gaussian_influence_vector!(out::AbstractVector, tel, dm, actuator_index::Int, fields::Tuple, p::AbstractVector)
    T = eltype(p)
    n = tel.params.resolution
    length(out) == n * n ||
        throw(DimensionMismatchError("Gaussian influence output length must match telescope resolution"))

    dm_topology = topology(dm)
    coords = actuator_coordinates(dm_topology)
    mis = dm.params.misregistration

    sigma = if :mechanical_coupling in fields
        coupling = _field_value(fields, p, :mechanical_coupling, mechanical_coupling(dm))
        _influence_width_from_mechanical_coupling(topology_axis_count(dm_topology), coupling)
    else
        _field_value(fields, p, :influence_width, influence_width(dm))
    end

    x0 = _field_value(fields, p, :actuator_x, coords[1, actuator_index])
    y0 = _field_value(fields, p, :actuator_y, coords[2, actuator_index])
    shift_x = _field_value(fields, p, :shift_x, mis.shift_x)
    shift_y = _field_value(fields, p, :shift_y, mis.shift_y)
    rot_deg = _field_value(fields, p, :rotation_deg, rotation_deg(mis))
    radial_scaling = _field_value(fields, p, :radial_scaling, mis.radial_scaling)
    tangential_scaling = _field_value(fields, p, :tangential_scaling, mis.tangential_scaling)

    x_m, y_m = _apply_ad_misregistration(
        x0,
        y0,
        shift_x,
        shift_y,
        rot_deg,
        convert(T, anamorphosis_angle_deg(mis)),
        tangential_scaling,
        radial_scaling,
    )

    cx = convert(T, (n + 1) / 2)
    cy = cx
    scale = convert(T, n / 2)
    two_sigma2 = convert(T, 2) * sigma^2
    @inbounds for j in 1:n
        y = (convert(T, j) - cy) / scale
        for i in 1:n
            x = (convert(T, i) - cx) / scale
            r2 = (x - x_m)^2 + (y - y_m)^2
            idx = (j - 1) * n + i
            out[idx] = ifelse(tel.state.pupil[i, j], exp(-r2 / two_sigma2), zero(T))
        end
    end
    return out
end

@inline function _influence_width_from_mechanical_coupling(n_act::Integer, coupling)
    T = typeof(coupling)
    pitch = convert(T, 2 / (n_act - 1))
    return pitch / sqrt(-convert(T, 2) * log(coupling))
end

function _apply_ad_misregistration(x, y, shift_x, shift_y, rotation_deg, anamorphosis_angle,
    tangential_scaling, radial_scaling)
    T = typeof(x + y + shift_x + shift_y + rotation_deg + anamorphosis_angle +
               tangential_scaling + radial_scaling)
    theta = convert(T, pi / 180) * anamorphosis_angle
    sin_theta, cos_theta = sincos(theta)
    phi = convert(T, pi / 180) * rotation_deg
    sin_phi, cos_phi = sincos(phi)

    b11 = tangential_scaling * cos_theta * cos_theta +
          radial_scaling * sin_theta * sin_theta
    b12 = (tangential_scaling - radial_scaling) * sin_theta * cos_theta
    b21 = b12
    b22 = tangential_scaling * sin_theta * sin_theta +
          radial_scaling * cos_theta * cos_theta

    m11 = cos_phi * b11 - sin_phi * b21
    m12 = cos_phi * b12 - sin_phi * b22
    m21 = sin_phi * b11 + cos_phi * b21
    m22 = sin_phi * b12 + cos_phi * b22
    xr = muladd(m12, y, m11 * x)
    yr = muladd(m22, y, m21 * x)
    return xr - shift_x, yr - shift_y
end

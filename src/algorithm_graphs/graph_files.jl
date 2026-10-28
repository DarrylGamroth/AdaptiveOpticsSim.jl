const ALGORITHM_GRAPH_SCHEMA_VERSION = 1

@inline function _file_error(context::AbstractString, message::AbstractString)
    return AlgorithmGraphError("$context: $message")
end

function _require_table_keys(
    table::AbstractDict,
    allowed::Tuple,
    required::Tuple,
    context::AbstractString,
)
    for key in keys(table)
        key in allowed || throw(_file_error(context, "unknown field '$key'"))
    end
    for key in required
        haskey(table, key) || throw(_file_error(context, "missing field '$key'"))
    end
    return nothing
end

function _require_named_fields(
    values::NamedTuple,
    allowed::Tuple,
    required::Tuple,
    context::AbstractString,
)
    for key in keys(values)
        key in allowed || throw(_file_error(context, "unknown field '$key'"))
    end
    for key in required
        hasproperty(values, key) || throw(_file_error(
            context,
            "missing field '$key'",
        ))
    end
    return nothing
end

function _file_string(value, context::AbstractString)
    value isa AbstractString || throw(_file_error(
        context,
        "expected a string, got $(typeof(value))",
    ))
    return String(value)
end

function _file_integer(value, context::AbstractString)
    value isa Integer && !(value isa Bool) || throw(_file_error(
        context,
        "expected an integer, got $(typeof(value))",
    ))
    return Int(value)
end

function _file_real(value, context::AbstractString)
    value isa Real && !(value isa Bool) || throw(_file_error(
        context,
        "expected a real number, got $(typeof(value))",
    ))
    isfinite(value) || throw(_file_error(context, "value must be finite"))
    return value
end

function _file_tuple(value, context::AbstractString)
    value isa Tuple || throw(_file_error(
        context,
        "expected an array, got $(typeof(value))",
    ))
    return value
end

function _file_real_tuple(value, context::AbstractString)
    values = _file_tuple(value, context)
    return ntuple(length(values)) do index
        _file_real(values[index], "$context[$index]")
    end
end

function _file_identifier_tuple(value, context::AbstractString)
    values = _file_tuple(value, context)
    return ntuple(length(values)) do index
        _file_identifier(values[index], "$context[$index]")
    end
end

function _file_bool(value, context::AbstractString)
    value isa Bool || throw(_file_error(
        context,
        "expected a boolean, got $(typeof(value))",
    ))
    return value
end

function _file_identifier(value, context::AbstractString)
    text = _file_string(value, context)
    Base.isidentifier(text) || throw(_file_error(
        context,
        "'$text' is not a valid identifier",
    ))
    return Symbol(text)
end

function _file_endpoint(value, context::AbstractString)
    text = _file_string(value, context)
    parts = split(text, '.'; keepempty=true)
    length(parts) == 2 || throw(_file_error(
        context,
        "endpoint '$text' must use node.port syntax",
    ))
    node = _file_identifier(parts[1], "$context node")
    port = _file_identifier(parts[2], "$context port")
    return AlgorithmPortReference(node, port)
end

_normalize_graph_value(value::Union{Bool,Integer,AbstractFloat,AbstractString}) =
    value

function _normalize_graph_value(values::AbstractVector)
    return Tuple(_normalize_graph_value(value) for value in values)
end


function _normalize_graph_value(table::AbstractDict)
    names = sort!(collect(String, keys(table)))
    symbols = Tuple(Symbol(name) for name in names)
    values = Tuple(_normalize_graph_value(table[name]) for name in names)
    return NamedTuple{symbols}(values)
end

function _normalize_graph_value(value)
    throw(AlgorithmGraphError(
        "unsupported graph value type $(typeof(value)) in config or props",
    ))
end

function _table_or_empty(table::AbstractDict, key::String, context::AbstractString)
    value = get(table, key, Dict{String,Union{}}())
    value isa AbstractDict || throw(_file_error(
        context,
        "field '$key' must be a table",
    ))
    return value
end

function _table_array_or_empty(
    table::AbstractDict,
    key::String,
    context::AbstractString,
)
    value = get(table, key, ())
    value isa AbstractVector || value isa Tuple || throw(_file_error(
        context,
        "field '$key' must be an array of tables",
    ))
    for (index, entry) in pairs(value)
        entry isa AbstractDict || throw(_file_error(
            context,
            "field '$key' entry $index must be a table",
        ))
    end
    return value
end

function _require_binding(
    bindings::NamedTuple,
    value,
    context::AbstractString,
)
    name = _file_identifier(value, context)
    hasproperty(bindings, name) || throw(_file_error(
        context,
        "no external binding named '$name' was supplied",
    ))
    return getproperty(bindings, name)
end

function _second_order_response_f32_node(name::Symbol, config::NamedTuple, props::NamedTuple)
    fields = (:natural_frequency_hz, :damping_ratio, :sample_period_s, :count, :schema)
    context = "node '$name' config"
    _require_named_fields(config, fields, fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    return second_order_response_node(name;
        natural_frequency_hz=_file_real(config.natural_frequency_hz, "$context.natural_frequency_hz"),
        damping_ratio=_file_real(config.damping_ratio, "$context.damping_ratio"),
        sample_period_s=_file_real(config.sample_period_s, "$context.sample_period_s"),
        count=_file_integer(config.count, "$context.count"),
        schema=_file_string(config.schema, "$context.schema"), T=Float32)
end

function _modal_opd_expansion_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :basis_schema,
        :coefficients_schema,
        :mode_count,
        :opd_schema,
        :pupil_columns,
        :pupil_rows,
        :pupil_support_schema,
    )
    _require_named_fields(
        config,
        config_fields,
        config_fields,
        "node '$name' config",
    )
    _require_named_fields(props, (), (), "node '$name' props")
    return modal_opd_expansion_node(
        name;
        pupil_rows=_file_integer(
            config.pupil_rows,
            "node '$name' config.pupil_rows",
        ),
        pupil_columns=_file_integer(
            config.pupil_columns,
            "node '$name' config.pupil_columns",
        ),
        mode_count=_file_integer(
            config.mode_count,
            "node '$name' config.mode_count",
        ),
        coefficients_schema=_file_string(
            config.coefficients_schema,
            "node '$name' config.coefficients_schema",
        ),
        opd_schema=_file_string(
            config.opd_schema,
            "node '$name' config.opd_schema",
        ),
        basis_schema=_file_string(
            config.basis_schema,
            "node '$name' config.basis_schema",
        ),
        pupil_support_schema=_file_string(
            config.pupil_support_schema,
            "node '$name' config.pupil_support_schema",
        ),
        T=Float32,
    )
end

function _deformable_mirror_surface_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :actuator_coordinates_schema,
        :actuator_count,
        :aperture_revision,
        :central_obstruction_ratio,
        :influence_functions_schema,
        :pdm_command_schema,
        :pupil_reflectivity,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    required_fields = (
        :actuator_coordinates_schema,
        :actuator_count,
        :influence_functions_schema,
        :pdm_command_schema,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    return deformable_mirror_surface_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count=_file_integer(
            config.actuator_count,
            "$context.actuator_count",
        ),
        pdm_command_schema=_file_string(
            config.pdm_command_schema,
            "$context.pdm_command_schema",
        ),
        surface_opd_schema=_file_string(
            config.surface_opd_schema,
            "$context.surface_opd_schema",
        ),
        actuator_coordinates_schema=_file_string(
            config.actuator_coordinates_schema,
            "$context.actuator_coordinates_schema",
        ),
        influence_functions_schema=_file_string(
            config.influence_functions_schema,
            "$context.influence_functions_schema",
        ),
        T=Float32,
    )
end

function _multilayer_atmosphere_opd_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :L0,
        :altitude,
        :aperture_revision,
        :atmosphere_opd_schema,
        :atmosphere_step,
        :central_obstruction_ratio,
        :fractional_cn2,
        :layer_ids,
        :pupil_reflectivity,
        :r0,
        :reference_wavelength_m,
        :resolution,
        :rng_seed,
        :source_separation_arcsec,
        :source_position_angle_deg,
        :telescope_diameter_m,
        :wind_direction_deg,
        :wind_speed,
    )
    required_fields = (
        :altitude,
        :atmosphere_opd_schema,
        :atmosphere_step,
        :fractional_cn2,
        :layer_ids,
        :r0,
        :reference_wavelength_m,
        :resolution,
        :rng_seed,
        :telescope_diameter_m,
        :wind_direction_deg,
        :wind_speed,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    L0 = hasproperty(config, :L0) ?
        _file_real(config.L0, "$context.L0") : 25.0
    return multilayer_atmosphere_opd_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        r0=_file_real(config.r0, "$context.r0"),
        reference_wavelength_m=_file_real(
            config.reference_wavelength_m,
            "$context.reference_wavelength_m",
        ),
        L0,
        source_separation_arcsec=hasproperty(config, :source_separation_arcsec) ?
            _file_real(config.source_separation_arcsec, "$context.source_separation_arcsec") : 0.0,
        source_position_angle_deg=hasproperty(config, :source_position_angle_deg) ?
            _file_real(config.source_position_angle_deg, "$context.source_position_angle_deg") : 0.0,
        fractional_cn2=_file_real_tuple(
            config.fractional_cn2,
            "$context.fractional_cn2",
        ),
        wind_speed=_file_real_tuple(
            config.wind_speed,
            "$context.wind_speed",
        ),
        wind_direction_deg=_file_real_tuple(
            config.wind_direction_deg,
            "$context.wind_direction_deg",
        ),
        altitude=_file_real_tuple(config.altitude, "$context.altitude"),
        layer_ids=_file_identifier_tuple(
            config.layer_ids,
            "$context.layer_ids",
        ),
        atmosphere_step=_file_real(
            config.atmosphere_step,
            "$context.atmosphere_step",
        ),
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        atmosphere_opd_schema=_file_string(
            config.atmosphere_opd_schema,
            "$context.atmosphere_opd_schema",
        ),
        T=Float32,
    )
end

function _spatial_ao_residual_f32_node(name::Symbol, config::NamedTuple,
        props::NamedTuple)
    fields = (:control_cutoff_cycles_per_m, :residual_amplitude, :residual_opd_schema)
    context = "node '$name' config"
    for field in fields
        hasproperty(config, field) || throw(_file_error(context, "missing field '$field'"))
    end
    hasproperty(config, :atmosphere_opd_schema) && throw(_file_error(context,
        "spatial AO residual uses residual_opd_schema, not atmosphere_opd_schema"))
    rejection = SpatialAORejection(
        control_cutoff_cycles_per_m=_file_real(config.control_cutoff_cycles_per_m,
            "$context.control_cutoff_cycles_per_m"),
        residual_amplitude=_file_real(config.residual_amplitude, "$context.residual_amplitude"))
    # Reuse the finite atmosphere's strict cold parser and defaults; only the
    # rejection and scientific output name differ from the raw declaration.
    atmosphere_config = merge(Base.structdiff(config, NamedTuple{fields}),
        (atmosphere_opd_schema=config.residual_opd_schema,))
    atmosphere = _multilayer_atmosphere_opd_f32_node(name, atmosphere_config, props)
    return algorithm_node(name, SpatialAOResidualNode{Float32},
        SpatialAOResidualNodeConfig(atmosphere.config, rejection); props)
end

function _gaussian_deformable_mirror_surface_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :actuator_coordinates_schema,
        :actuator_count,
        :aperture_revision,
        :central_obstruction_ratio,
        :influence_width,
        :pdm_command_schema,
        :pupil_reflectivity,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    required_fields = (
        :actuator_coordinates_schema,
        :actuator_count,
        :influence_width,
        :pdm_command_schema,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    return gaussian_deformable_mirror_surface_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count=_file_integer(
            config.actuator_count,
            "$context.actuator_count",
        ),
        influence_width=_file_real(
            config.influence_width,
            "$context.influence_width",
        ),
        pdm_command_schema=_file_string(
            config.pdm_command_schema,
            "$context.pdm_command_schema",
        ),
        surface_opd_schema=_file_string(
            config.surface_opd_schema,
            "$context.surface_opd_schema",
        ),
        actuator_coordinates_schema=_file_string(
            config.actuator_coordinates_schema,
            "$context.actuator_coordinates_schema",
        ),
        T=Float32,
    )
end

function _grid_gaussian_deformable_mirror_surface_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :actuator_axis_count,
        :actuator_count,
        :actuator_grid_indices_schema,
        :actuator_pitch,
        :aperture_revision,
        :central_obstruction_ratio,
        :influence_width,
        :pdm_command_schema,
        :pupil_reflectivity,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    required_fields = (
        :actuator_axis_count,
        :actuator_count,
        :actuator_grid_indices_schema,
        :actuator_pitch,
        :influence_width,
        :pdm_command_schema,
        :resolution,
        :surface_opd_schema,
        :telescope_diameter_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    return grid_gaussian_deformable_mirror_surface_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count=_file_integer(
            config.actuator_count,
            "$context.actuator_count",
        ),
        actuator_axis_count=_file_integer(
            config.actuator_axis_count,
            "$context.actuator_axis_count",
        ),
        actuator_pitch=_file_real(
            config.actuator_pitch,
            "$context.actuator_pitch",
        ),
        influence_width=_file_real(
            config.influence_width,
            "$context.influence_width",
        ),
        pdm_command_schema=_file_string(
            config.pdm_command_schema,
            "$context.pdm_command_schema",
        ),
        surface_opd_schema=_file_string(
            config.surface_opd_schema,
            "$context.surface_opd_schema",
        ),
        actuator_grid_indices_schema=_file_string(
            config.actuator_grid_indices_schema,
            "$context.actuator_grid_indices_schema",
        ),
        T=Float32,
    )
end

function _pupil_opd_composition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :pupil_opd_schema,
        :resolution,
        :surface_opd_schema,
        :uncompensated_opd_schema,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, config_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    return pupil_opd_composition_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        uncompensated_opd_schema=_file_string(
            config.uncompensated_opd_schema,
            "$context.uncompensated_opd_schema",
        ),
        surface_opd_schema=_file_string(
            config.surface_opd_schema,
            "$context.surface_opd_schema",
        ),
        pupil_opd_schema=_file_string(
            config.pupil_opd_schema,
            "$context.pupil_opd_schema",
        ),
        T=Float32,
    )
end

function _shack_hartmann_rate_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :aperture_revision,
        :central_obstruction_ratio,
        :diffraction_padding,
        :half_pixel_shift,
        :n_lenslets,
        :n_pix_subap,
        :opd_schema,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :pupil_reflectivity,
        :resolution,
        :shannon_sampling,
        :source_band,
        :source_magnitude,
        :source_photon_irradiance_m2_s,
        :source_position_angle_deg,
        :source_separation_arcsec,
        :source_wavelength_m,
        :telescope_diameter_m,
        :threshold_convolution,
        :valid_subaperture_threshold,
    )
    required_fields = (
        :n_lenslets,
        :n_pix_subap,
        :opd_schema,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :resolution,
        :source_wavelength_m,
        :telescope_diameter_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")

    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    diffraction_padding = hasproperty(config, :diffraction_padding) ?
        _file_integer(
            config.diffraction_padding,
            "$context.diffraction_padding",
        ) : 2
    valid_subaperture_threshold = hasproperty(
        config,
        :valid_subaperture_threshold,
    ) ? _file_real(
        config.valid_subaperture_threshold,
        "$context.valid_subaperture_threshold",
    ) : 0.1
    threshold_convolution = hasproperty(config, :threshold_convolution) ?
        _file_real(
            config.threshold_convolution,
            "$context.threshold_convolution",
        ) : 0.05
    half_pixel_shift = hasproperty(config, :half_pixel_shift) ?
        _file_bool(config.half_pixel_shift, "$context.half_pixel_shift") :
        false
    shannon_sampling = hasproperty(config, :shannon_sampling) ?
        _file_bool(config.shannon_sampling, "$context.shannon_sampling") :
        true
    source_separation_arcsec = hasproperty(
        config,
        :source_separation_arcsec,
    ) ? _file_real(
        config.source_separation_arcsec,
        "$context.source_separation_arcsec",
    ) : 0.0
    source_position_angle_deg = hasproperty(
        config,
        :source_position_angle_deg,
    ) ? _file_real(
        config.source_position_angle_deg,
        "$context.source_position_angle_deg",
    ) : 0.0
    source_band = hasproperty(config, :source_band) ? _file_identifier(
        config.source_band,
        "$context.source_band",
    ) : :custom
    source_magnitude = hasproperty(config, :source_magnitude) ? _file_real(
        config.source_magnitude,
        "$context.source_magnitude",
    ) : 0.0
    source_photon_irradiance_m2_s = hasproperty(
        config,
        :source_photon_irradiance_m2_s,
    ) ? _file_real(
        config.source_photon_irradiance_m2_s,
        "$context.source_photon_irradiance_m2_s",
    ) : nothing

    return shack_hartmann_rate_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        n_lenslets=_file_integer(config.n_lenslets, "$context.n_lenslets"),
        n_pix_subap=_file_integer(
            config.n_pix_subap,
            "$context.n_pix_subap",
        ),
        diffraction_padding,
        pixel_scale_arcsec=_file_real(
            config.pixel_scale_arcsec,
            "$context.pixel_scale_arcsec",
        ),
        valid_subaperture_threshold,
        threshold_convolution,
        half_pixel_shift,
        shannon_sampling,
        source_band,
        source_magnitude,
        source_wavelength_m=_file_real(
            config.source_wavelength_m,
            "$context.source_wavelength_m",
        ),
        source_photon_irradiance_m2_s,
        source_separation_arcsec,
        source_position_angle_deg,
        opd_schema=_file_string(config.opd_schema, "$context.opd_schema"),
        photon_rate_schema=_file_string(
            config.photon_rate_schema,
            "$context.photon_rate_schema",
        ),
        T=Float32,
    )
end

function _pyramid_rate_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :aperture_revision,
        :binning,
        :central_obstruction_ratio,
        :diffraction_padding,
        :modulation,
        :modulation_points,
        :modulation_propagation_strategy,
        :n_pix_edge,
        :n_pix_separation,
        :opd_schema,
        :photon_rate_schema,
        :psf_centering,
        :pupil_reflectivity,
        :pupil_samples,
        :pupil_shift_x_pixels,
        :pupil_shift_y_pixels,
        :resolution,
        :source_band,
        :source_magnitude,
        :source_photon_irradiance_m2_s,
        :source_position_angle_deg,
        :source_separation_arcsec,
        :source_wavelength_m,
        :telescope_diameter_m,
    )
    required_fields = (
        :opd_schema,
        :photon_rate_schema,
        :pupil_samples,
        :resolution,
        :source_wavelength_m,
        :telescope_diameter_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")

    central_obstruction_ratio = hasproperty(
        config,
        :central_obstruction_ratio,
    ) ? _file_real(
        config.central_obstruction_ratio,
        "$context.central_obstruction_ratio",
    ) : 0.0
    pupil_reflectivity = hasproperty(config, :pupil_reflectivity) ?
        _file_real(config.pupil_reflectivity, "$context.pupil_reflectivity") :
        1.0
    aperture_revision = hasproperty(config, :aperture_revision) ?
        _file_integer(config.aperture_revision, "$context.aperture_revision") :
        0
    modulation = hasproperty(config, :modulation) ?
        _file_real(config.modulation, "$context.modulation") : 2.0
    modulation_points = hasproperty(config, :modulation_points) ?
        _file_integer(
            config.modulation_points,
            "$context.modulation_points",
        ) : nothing
    modulation_propagation_strategy = if hasproperty(
        config,
        :modulation_propagation_strategy,
    )
        strategy_name = _file_identifier(
            config.modulation_propagation_strategy,
            "$context.modulation_propagation_strategy",
        )
        if strategy_name === :pupil_tilt
            PyramidPupilTiltStrategy()
        else
            throw(AlgorithmGraphError(
                "$context.modulation_propagation_strategy must be " *
                "'pupil_tilt'",
            ))
        end
    else
        PyramidPupilTiltStrategy()
    end
    diffraction_padding = hasproperty(config, :diffraction_padding) ?
        _file_integer(
            config.diffraction_padding,
            "$context.diffraction_padding",
        ) : 2
    psf_centering = hasproperty(config, :psf_centering) ?
        _file_bool(config.psf_centering, "$context.psf_centering") : true
    n_pix_separation = hasproperty(config, :n_pix_separation) ?
        _file_integer(
            config.n_pix_separation,
            "$context.n_pix_separation",
        ) : nothing
    n_pix_edge = hasproperty(config, :n_pix_edge) ?
        _file_integer(config.n_pix_edge, "$context.n_pix_edge") : nothing
    binning = hasproperty(config, :binning) ?
        _file_integer(config.binning, "$context.binning") : 1
    pupil_shifts = map((:pupil_shift_x_pixels, :pupil_shift_y_pixels)) do field
        hasproperty(config, field) || return 0
        value = getproperty(config, field)
        field_context = "$context.$field"
        return value isa Real ? _file_real(value, field_context) :
            _file_real_tuple(value, field_context)
    end
    source_band = hasproperty(config, :source_band) ? _file_identifier(
        config.source_band,
        "$context.source_band",
    ) : :custom
    source_magnitude = hasproperty(config, :source_magnitude) ? _file_real(
        config.source_magnitude,
        "$context.source_magnitude",
    ) : 0.0
    source_photon_irradiance_m2_s = hasproperty(
        config,
        :source_photon_irradiance_m2_s,
    ) ? _file_real(
        config.source_photon_irradiance_m2_s,
        "$context.source_photon_irradiance_m2_s",
    ) : nothing
    source_separation_arcsec = hasproperty(
        config,
        :source_separation_arcsec,
    ) ? _file_real(
        config.source_separation_arcsec,
        "$context.source_separation_arcsec",
    ) : 0.0
    source_position_angle_deg = hasproperty(
        config,
        :source_position_angle_deg,
    ) ? _file_real(
        config.source_position_angle_deg,
        "$context.source_position_angle_deg",
    ) : 0.0

    return pyramid_rate_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        pupil_samples=_file_integer(
            config.pupil_samples,
            "$context.pupil_samples",
        ),
        modulation,
        modulation_points,
        modulation_propagation_strategy,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge,
        binning,
        pupil_shift_x_pixels=pupil_shifts[1],
        pupil_shift_y_pixels=pupil_shifts[2],
        source_band,
        source_magnitude,
        source_wavelength_m=_file_real(
            config.source_wavelength_m,
            "$context.source_wavelength_m",
        ),
        source_photon_irradiance_m2_s,
        source_separation_arcsec,
        source_position_angle_deg,
        opd_schema=_file_string(config.opd_schema, "$context.opd_schema"),
        photon_rate_schema=_file_string(
            config.photon_rate_schema,
            "$context.photon_rate_schema",
        ),
        T=Float32,
    )
end

function _ccd_detector_acquisition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :binning,
        :columns,
        :exposure_duration_s,
        :frame_schema,
        :photon_noise,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :quantum_efficiency,
        :readout_noise,
        :readout_noise_e,
        :rng_seed,
        :rows,
        :wavelength_m,
    )
    required_fields = (
        :columns,
        :exposure_duration_s,
        :frame_schema,
        :photon_noise,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :quantum_efficiency,
        :readout_noise,
        :readout_noise_e,
        :rng_seed,
        :rows,
        :wavelength_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    binning = hasproperty(config, :binning) ? _file_integer(
        config.binning,
        "$context.binning",
    ) : 1
    return ccd_detector_acquisition_node(
        name;
        rows=_file_integer(config.rows, "$context.rows"),
        columns=_file_integer(config.columns, "$context.columns"),
        binning,
        pixel_scale_arcsec=_file_real(
            config.pixel_scale_arcsec,
            "$context.pixel_scale_arcsec",
        ),
        wavelength_m=_file_real(
            config.wavelength_m,
            "$context.wavelength_m",
        ),
        exposure_duration_s=_file_real(
            config.exposure_duration_s,
            "$context.exposure_duration_s",
        ),
        quantum_efficiency=_file_real(
            config.quantum_efficiency,
            "$context.quantum_efficiency",
        ),
        photon_noise=_file_bool(
            config.photon_noise,
            "$context.photon_noise",
        ),
        readout_noise=_file_bool(
            config.readout_noise,
            "$context.readout_noise",
        ),
        readout_noise_e=_file_real(
            config.readout_noise_e,
            "$context.readout_noise_e",
        ),
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        photon_rate_schema=_file_string(
            config.photon_rate_schema,
            "$context.photon_rate_schema",
        ),
        frame_schema=_file_string(
            config.frame_schema,
            "$context.frame_schema",
        ),
        T=Float32,
    )
end

function _cmos_detector_acquisition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :binning,
        :bits,
        :column_readout_noise_e,
        :columns,
        :dark_current_e_per_pixel_s,
        :exposure_duration_s,
        :frame_schema,
        :full_well_e,
        :gain,
        :photon_noise,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :quantum_efficiency,
        :readout_noise,
        :readout_noise_e,
        :rng_seed,
        :row_readout_noise_e,
        :rows,
        :wavelength_m,
    )
    required_fields = (
        :columns,
        :exposure_duration_s,
        :frame_schema,
        :photon_rate_schema,
        :pixel_scale_arcsec,
        :quantum_efficiency,
        :rng_seed,
        :rows,
        :wavelength_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")

    binning = hasproperty(config, :binning) ? _file_integer(
        config.binning,
        "$context.binning",
    ) : 1
    gain = hasproperty(config, :gain) ?
        _file_real(config.gain, "$context.gain") : 1.0
    dark_current = hasproperty(config, :dark_current_e_per_pixel_s) ?
        _file_real(
            config.dark_current_e_per_pixel_s,
            "$context.dark_current_e_per_pixel_s",
        ) : 0.0
    bits = hasproperty(config, :bits) ?
        _file_integer(config.bits, "$context.bits") : nothing
    full_well = hasproperty(config, :full_well_e) ?
        _file_real(config.full_well_e, "$context.full_well_e") : nothing
    photon_noise = hasproperty(config, :photon_noise) ?
        _file_bool(config.photon_noise, "$context.photon_noise") : true
    readout_noise = hasproperty(config, :readout_noise) ?
        _file_bool(config.readout_noise, "$context.readout_noise") : false
    readout_noise_e = hasproperty(config, :readout_noise_e) ?
        _file_real(config.readout_noise_e, "$context.readout_noise_e") : 0.0
    column_readout_noise = hasproperty(config, :column_readout_noise_e) ?
        _file_real(
            config.column_readout_noise_e,
            "$context.column_readout_noise_e",
        ) : 0.0
    row_readout_noise = hasproperty(config, :row_readout_noise_e) ?
        _file_real(
            config.row_readout_noise_e,
            "$context.row_readout_noise_e",
        ) : 0.0

    return cmos_detector_acquisition_node(
        name;
        rows=_file_integer(config.rows, "$context.rows"),
        columns=_file_integer(config.columns, "$context.columns"),
        binning,
        pixel_scale_arcsec=_file_real(
            config.pixel_scale_arcsec,
            "$context.pixel_scale_arcsec",
        ),
        wavelength_m=_file_real(
            config.wavelength_m,
            "$context.wavelength_m",
        ),
        exposure_duration_s=_file_real(
            config.exposure_duration_s,
            "$context.exposure_duration_s",
        ),
        quantum_efficiency=_file_real(
            config.quantum_efficiency,
            "$context.quantum_efficiency",
        ),
        gain,
        dark_current_e_per_pixel_s=dark_current,
        bits,
        full_well_e=full_well,
        photon_noise,
        readout_noise,
        readout_noise_e,
        column_readout_noise_e=column_readout_noise,
        row_readout_noise_e=row_readout_noise,
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        photon_rate_schema=_file_string(
            config.photon_rate_schema,
            "$context.photon_rate_schema",
        ),
        frame_schema=_file_string(
            config.frame_schema,
            "$context.frame_schema",
        ),
        T=Float32,
    )
end

function _emccd_detector_acquisition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :binning,
        :bits,
        :clock_induced_charge_e_per_pixel_frame,
        :columns,
        :dark_current_e_per_pixel_s,
        :em_gain_max,
        :em_gain_min,
        :excess_noise_factor,
        :exposure_duration_s,
        :frame_schema,
        :full_well_e,
        :gain,
        :normalized_pupil_sampling,
        :photon_noise,
        :photon_rate_schema,
        :quantum_efficiency,
        :readout_noise,
        :readout_noise_e,
        :register_full_well_e,
        :rng_seed,
        :rows,
        :wavelength_m,
    )
    required_fields = (
        :columns,
        :exposure_duration_s,
        :frame_schema,
        :normalized_pupil_sampling,
        :photon_rate_schema,
        :quantum_efficiency,
        :rng_seed,
        :rows,
        :wavelength_m,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")

    binning = hasproperty(config, :binning) ? _file_integer(
        config.binning,
        "$context.binning",
    ) : 1
    gain = hasproperty(config, :gain) ?
        _file_real(config.gain, "$context.gain") : 1.0
    dark_current = hasproperty(config, :dark_current_e_per_pixel_s) ?
        _file_real(
            config.dark_current_e_per_pixel_s,
            "$context.dark_current_e_per_pixel_s",
        ) : 0.0
    bits = hasproperty(config, :bits) ?
        _file_integer(config.bits, "$context.bits") : nothing
    full_well = hasproperty(config, :full_well_e) ?
        _file_real(config.full_well_e, "$context.full_well_e") : nothing
    photon_noise = hasproperty(config, :photon_noise) ?
        _file_bool(config.photon_noise, "$context.photon_noise") : true
    readout_noise = hasproperty(config, :readout_noise) ?
        _file_bool(config.readout_noise, "$context.readout_noise") : false
    readout_noise_e = hasproperty(config, :readout_noise_e) ?
        _file_real(config.readout_noise_e, "$context.readout_noise_e") : 0.0
    excess_noise_factor = hasproperty(config, :excess_noise_factor) ?
        _file_real(
            config.excess_noise_factor,
            "$context.excess_noise_factor",
        ) : 1.0
    cic = hasproperty(
        config,
        :clock_induced_charge_e_per_pixel_frame,
    ) ? _file_real(
        config.clock_induced_charge_e_per_pixel_frame,
        "$context.clock_induced_charge_e_per_pixel_frame",
    ) : 0.0
    register_full_well = hasproperty(config, :register_full_well_e) ?
        _file_real(
            config.register_full_well_e,
            "$context.register_full_well_e",
        ) : nothing
    em_gain_min = hasproperty(config, :em_gain_min) ?
        _file_real(config.em_gain_min, "$context.em_gain_min") : 1.0
    em_gain_max = hasproperty(config, :em_gain_max) ?
        _file_real(config.em_gain_max, "$context.em_gain_max") : Inf

    return emccd_detector_acquisition_node(
        name;
        rows=_file_integer(config.rows, "$context.rows"),
        columns=_file_integer(config.columns, "$context.columns"),
        binning,
        normalized_pupil_sampling=_file_real(
            config.normalized_pupil_sampling,
            "$context.normalized_pupil_sampling",
        ),
        wavelength_m=_file_real(
            config.wavelength_m,
            "$context.wavelength_m",
        ),
        exposure_duration_s=_file_real(
            config.exposure_duration_s,
            "$context.exposure_duration_s",
        ),
        quantum_efficiency=_file_real(
            config.quantum_efficiency,
            "$context.quantum_efficiency",
        ),
        gain,
        dark_current_e_per_pixel_s=dark_current,
        bits,
        full_well_e=full_well,
        photon_noise,
        readout_noise,
        readout_noise_e,
        excess_noise_factor,
        clock_induced_charge_e_per_pixel_frame=cic,
        register_full_well_e=register_full_well,
        em_gain_min,
        em_gain_max,
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        photon_rate_schema=_file_string(
            config.photon_rate_schema,
            "$context.photon_rate_schema",
        ),
        frame_schema=_file_string(
            config.frame_schema,
            "$context.frame_schema",
        ),
        T=Float32,
    )
end

_file_persistence_model(value, context::AbstractString) =
    throw(_file_error(context, "expected a mapping"))

function _file_persistence_model(value::NamedTuple, context::AbstractString)
    _require_named_fields(value, (:mode, :coupling, :decay), (:mode,), context)
    mode = _file_string(value.mode, "$context.mode")
    if mode == "none"
        _require_named_fields(value, (:mode,), (:mode,), context)
        return NullPersistence()
    elseif mode == "exponential"
        fields = (:mode, :coupling, :decay)
        _require_named_fields(value, fields, fields, context)
        coupling = _file_real(value.coupling, "$context.coupling")
        decay = _file_real(value.decay, "$context.decay")
        0 <= coupling <= 1 || throw(_file_error(context, "coupling must lie in [0, 1]"))
        0 <= decay <= 1 || throw(_file_error(context, "decay must lie in [0, 1]"))
        return ExponentialPersistence(coupling, decay)
    end
    throw(_file_error(context, "unsupported persistence mode '$mode'"))
end

_file_nonlinearity_model(value, context::AbstractString) =
    throw(_file_error(context, "expected a mapping"))

function _file_nonlinearity_model(value::NamedTuple, context::AbstractString)
    _require_named_fields(value, (:mode, :coefficient), (:mode,), context)
    mode = _file_string(value.mode, "$context.mode")
    if mode == "none"
        _require_named_fields(value, (:mode,), (:mode,), context)
        return NullFrameNonlinearity()
    elseif mode == "saturating"
        fields = (:mode, :coefficient)
        _require_named_fields(value, fields, fields, context)
        coefficient = _file_real(value.coefficient, "$context.coefficient")
        coefficient >= 0 || throw(_file_error(context, "coefficient must be nonnegative"))
        return SaturatingFrameNonlinearity(coefficient)
    end
    throw(_file_error(context, "unsupported nonlinearity mode '$mode'"))
end

function _ingaas_detector_acquisition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    required_fields = (:rows, :columns, :pixel_sampling_m, :spectral_channel,
        :exposure_duration_s, :quantum_efficiency, :rng_seed,
        :photon_rate_schema, :frame_schema)
    config_fields = (required_fields..., :gain, :dark_current_e_per_pixel_s,
        :glow_e_per_pixel_s, :bits, :full_well_e, :photon_noise, :readout_noise,
        :readout_noise_e, :persistence_model, :nonlinearity_model)
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    sampling = _file_real_tuple(config.pixel_sampling_m, "$context.pixel_sampling_m")
    length(sampling) == 2 || throw(_file_error(
        "$context.pixel_sampling_m", "expected row and column pitches"))
    return ingaas_detector_acquisition_node(name;
        rows=_file_integer(config.rows, "$context.rows"),
        columns=_file_integer(config.columns, "$context.columns"),
        pixel_sampling_m=sampling,
        spectral_channel=_file_identifier(config.spectral_channel, "$context.spectral_channel"),
        exposure_duration_s=_file_real(config.exposure_duration_s, "$context.exposure_duration_s"),
        quantum_efficiency=_file_real(config.quantum_efficiency, "$context.quantum_efficiency"),
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        photon_rate_schema=_file_string(config.photon_rate_schema, "$context.photon_rate_schema"),
        frame_schema=_file_string(config.frame_schema, "$context.frame_schema"),
        gain=_file_real(get(config, :gain, 1), "$context.gain"),
        dark_current_e_per_pixel_s=_file_real(get(config, :dark_current_e_per_pixel_s, 0),
            "$context.dark_current_e_per_pixel_s"),
        glow_e_per_pixel_s=_file_real(get(config, :glow_e_per_pixel_s, 0),
            "$context.glow_e_per_pixel_s"),
        bits=hasproperty(config, :bits) ? _file_integer(config.bits, "$context.bits") : nothing,
        full_well_e=hasproperty(config, :full_well_e) ?
            _file_real(config.full_well_e, "$context.full_well_e") : nothing,
        photon_noise=_file_bool(get(config, :photon_noise, true), "$context.photon_noise"),
        readout_noise=_file_bool(get(config, :readout_noise, false), "$context.readout_noise"),
        readout_noise_e=_file_real(get(config, :readout_noise_e, 0), "$context.readout_noise_e"),
        persistence_model=hasproperty(config, :persistence_model) ?
            _file_persistence_model(config.persistence_model, "$context.persistence_model") :
            NullPersistence(),
        nonlinearity_model=hasproperty(config, :nonlinearity_model) ?
            _file_nonlinearity_model(config.nonlinearity_model, "$context.nonlinearity_model") :
            NullFrameNonlinearity(),
        T=Float32)
end

function _file_hgcdte_avalanche_multiplication(value, context::AbstractString)
    model = _file_string(value, context)
    model == "conditional_gamma" && return ConditionalGammaAvalancheMultiplication()
    model == "clipped_gaussian_approximation" &&
        return ClippedGaussianAvalancheMultiplicationApproximation()
    throw(_file_error(context, "unsupported avalanche multiplication model '$model'"))
end

function _hgcdte_avalanche_detector_acquisition_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    required_fields = (:rows, :columns, :pixel_sampling_m, :spectral_channel,
        :exposure_duration_s, :quantum_efficiency, :rng_seed,
        :photon_rate_schema, :frame_schema)
    config_fields = (required_fields..., :gain, :dark_current_e_per_pixel_s,
        :glow_e_per_pixel_s, :bits, :full_well_e, :photon_noise, :readout_noise,
        :readout_noise_e, :persistence_model, :nonlinearity_model,
        :avalanche_gain, :excess_noise_factor, :multiplication_model)
    context = "node '$name' config"
    _require_named_fields(config, config_fields, required_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    sampling = _file_real_tuple(config.pixel_sampling_m, "$context.pixel_sampling_m")
    length(sampling) == 2 || throw(_file_error(
        "$context.pixel_sampling_m", "expected row and column pitches"))
    return hgcdte_avalanche_detector_acquisition_node(name;
        rows=_file_integer(config.rows, "$context.rows"),
        columns=_file_integer(config.columns, "$context.columns"),
        pixel_sampling_m=sampling,
        spectral_channel=_file_identifier(config.spectral_channel, "$context.spectral_channel"),
        exposure_duration_s=_file_real(config.exposure_duration_s, "$context.exposure_duration_s"),
        quantum_efficiency=_file_real(config.quantum_efficiency, "$context.quantum_efficiency"),
        rng_seed=_file_integer(config.rng_seed, "$context.rng_seed"),
        photon_rate_schema=_file_string(config.photon_rate_schema, "$context.photon_rate_schema"),
        frame_schema=_file_string(config.frame_schema, "$context.frame_schema"),
        gain=_file_real(get(config, :gain, 1), "$context.gain"),
        dark_current_e_per_pixel_s=_file_real(get(config, :dark_current_e_per_pixel_s, 0),
            "$context.dark_current_e_per_pixel_s"),
        glow_e_per_pixel_s=_file_real(get(config, :glow_e_per_pixel_s, 0),
            "$context.glow_e_per_pixel_s"),
        bits=hasproperty(config, :bits) ? _file_integer(config.bits, "$context.bits") : nothing,
        full_well_e=hasproperty(config, :full_well_e) ?
            _file_real(config.full_well_e, "$context.full_well_e") : nothing,
        photon_noise=_file_bool(get(config, :photon_noise, true), "$context.photon_noise"),
        readout_noise=_file_bool(get(config, :readout_noise, false), "$context.readout_noise"),
        readout_noise_e=_file_real(get(config, :readout_noise_e, 0), "$context.readout_noise_e"),
        persistence_model=hasproperty(config, :persistence_model) ?
            _file_persistence_model(config.persistence_model, "$context.persistence_model") :
            NullPersistence(),
        nonlinearity_model=hasproperty(config, :nonlinearity_model) ?
            _file_nonlinearity_model(config.nonlinearity_model, "$context.nonlinearity_model") :
            NullFrameNonlinearity(),
        avalanche_gain=_file_real(get(config, :avalanche_gain, 1), "$context.avalanche_gain"),
        excess_noise_factor=_file_real(get(config, :excess_noise_factor, 1),
            "$context.excess_noise_factor"),
        multiplication_model=_file_hgcdte_avalanche_multiplication(
            get(config, :multiplication_model, "clipped_gaussian_approximation"),
            "$context.multiplication_model"),
        T=Float32)
end

_file_matrix(value::AbstractMatrix, context::AbstractString) = value
_file_matrix(value, context::AbstractString) =
    throw(_file_error(context, "external binding must be a matrix"))

"""
    graph_node_bindings(factory, config::NamedTuple)

Return the tuple of external binding names (`Symbol`s) required by a cold node
factory's construction configuration. The default is `()`. A factory requiring
calibration arrays must specialize this method so applications can load only
those arrays before invoking `(name, config, props, bindings)`.

This declaration performs no I/O and owns no arrays. It does not describe frame
inputs, delayed initial values or sparse parameters; those have explicit graph
entries. The native constructor still validates and snapshots the calibration.
"""
graph_node_bindings(factory, config::NamedTuple) = ()

"""
    graph_node_configuration(factory, config::NamedTuple, binding_names::NamedTuple=NamedTuple())

Validate a node's explicit `bindings` mapping and return the factory configuration
with each binding role mapped to its external name. The names must match
[`graph_node_bindings`](@ref). Configuration keys and binding roles must not
collide. This cold operation loads no arrays and performs no I/O.
"""
function graph_node_configuration(factory, config::NamedTuple,
    binding_names::NamedTuple=NamedTuple())
    declared = map(keys(binding_names)) do role
        _file_identifier(String(role), "node binding role")
        hasproperty(config, role) && throw(_file_error(
            "node binding '$role'", "binding role also appears in config"))
        _file_identifier(getproperty(binding_names, role), "node binding '$role'")
    end
    merged = merge(config, binding_names)
    required = _file_binding_names(graph_node_bindings(factory, merged))
    Set(declared) == Set(required) || throw(AlgorithmGraphError(
        "node bindings must declare exactly the construction calibration names required by the factory"))
    return merged
end

_file_binding_names(names::Tuple{Vararg{Symbol}}) = names
_file_binding_names(names) = throw(AlgorithmGraphError(
    "graph_node_bindings must return a tuple of binding-name Symbols"))

function _pixel_mtf_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
    bindings::NamedTuple,
)
    context = "node '$name' config"
    fields = (:input_shape, :padding, :input_schema, :output_schema, :transfer)
    _require_named_fields(config, fields, fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    shape = _file_tuple(config.input_shape, "$context.input_shape")
    length(shape) == 2 || throw(_file_error("$context.input_shape", "expected two dimensions"))
    input_shape = ntuple(i -> _file_integer(shape[i], "$context.input_shape[$i]"), 2)
    transfer = _file_matrix(
        _require_binding(bindings, config.transfer, "$context.transfer"),
        "$context.transfer")
    return pixel_mtf_node(name; input_shape,
        padding=_file_integer(config.padding, "$context.padding"),
        input_schema=_file_string(config.input_schema, "$context.input_schema"),
        output_schema=_file_string(config.output_schema, "$context.output_schema"),
        transfer, T=Float32)
end

function graph_node_bindings(::typeof(_pixel_mtf_f32_node), config::NamedTuple)
    context = "pixel MTF config"
    fields = (:input_shape, :padding, :input_schema, :output_schema, :transfer)
    _require_named_fields(config, fields, fields, context)
    return (_file_identifier(config.transfer, "$context.transfer"),)
end

function _detector_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple;
    detector::Symbol,
)
    factories = (
        ccd=_ccd_detector_acquisition_f32_node,
        cmos=_cmos_detector_acquisition_f32_node,
        emccd=_emccd_detector_acquisition_f32_node,
        ingaas=_ingaas_detector_acquisition_f32_node,
        hgcdte_linear_avalanche_array=_hgcdte_avalanche_detector_acquisition_f32_node,
    )
    hasproperty(factories, detector) || throw(_file_error(
        "node '$name' detector",
        "unsupported detector '$detector'; supported detectors are ccd, cmos, emccd, ingaas, hgcdte_linear_avalanche_array",
    ))
    if hasproperty(config, :shutter)
        shutter = _file_string(config.shutter, "node '$name' config.shutter")
        shutter == "global" || throw(_file_error(
            "node '$name' config.shutter",
            "unsupported shutter '$shutter'; this graph node supports only global acquisition",
        ))
    end
    acquisition_config = Base.structdiff(config, NamedTuple{(:shutter,)})
    return getproperty(factories, detector)(name, acquisition_config, props)
end

"""
    builtin_graph_node_types()

Return the immutable node-type map accepted by the maintained declarative graph loader.
Built-in optional factories report their missing dependency when constructed.
Applications may merge additional cold factories into this named tuple
and pass the result explicitly to [`load_algorithm_graph`](@ref). No mutable
global registry or dynamic Julia evaluation is used.
"""
function builtin_graph_node_types()
    return (
        proper_optical_train_f64=proper_optical_train_node_factory(; T=Float64),
        proper_optical_train_f32=proper_optical_train_node_factory(; T=Float32),
        detector=_detector_node,
        pixel_mtf=_pixel_mtf_f32_node,
        deformable_mirror_surface_f32=
            _deformable_mirror_surface_f32_node,
        gaussian_deformable_mirror_surface_f32=
            _gaussian_deformable_mirror_surface_f32_node,
        grid_gaussian_deformable_mirror_surface_f32=
            _grid_gaussian_deformable_mirror_surface_f32_node,
        modal_opd_expansion_f32=_modal_opd_expansion_f32_node,
        second_order_response_f32=_second_order_response_f32_node,
        multilayer_atmosphere_opd_f32=
            _multilayer_atmosphere_opd_f32_node,
        spatial_ao_residual_f32=_spatial_ao_residual_f32_node,
        pupil_opd_composition_f32=_pupil_opd_composition_f32_node,
        pyramid_rate_f32=_pyramid_rate_f32_node,
        shack_hartmann_rate_f32=_shack_hartmann_rate_f32_node,
    )
end

function _parse_file_node(
    entry::AbstractDict,
    node_types::NamedTuple,
    bindings::NamedTuple,
    index::Int,
)
    context = "nodes[$index]"
    _require_table_keys(
        entry,
        ("name", "type", "detector", "config", "props", "bindings"),
        ("name", "type", "config"),
        context,
    )
    name = _file_identifier(entry["name"], "$context.name")
    node_type = _file_identifier(entry["type"], "$context.type")
    hasproperty(node_types, node_type) || throw(_file_error(
        "$context.type",
        "unknown graph node type '$node_type'",
    ))
    if node_type === :detector
        haskey(entry, "detector") || throw(_file_error(
            context,
            "missing field 'detector'",
        ))
    elseif haskey(entry, "detector")
        throw(_file_error(context, "field 'detector' is only valid for type 'detector'"))
    end
    config = _normalize_graph_value(
        _table_or_empty(entry, "config", context),
    )
    props = _normalize_graph_value(
        _table_or_empty(entry, "props", context),
    )
    factory = getproperty(node_types, node_type)
    config = graph_node_configuration(factory, config,
        _normalize_graph_value(_table_or_empty(entry, "bindings", context)))
    definition = if node_type === :detector
        detector = _file_identifier(entry["detector"], "$context.detector")
        factory(name, config, props; detector)
    elseif applicable(factory, name, config, props, bindings)
        factory(name, config, props, bindings)
    else
        factory(name, config, props)
    end
    definition isa AlgorithmNodeDefinition || throw(_file_error(
        context,
        "node type '$node_type' did not return an AlgorithmNodeDefinition",
    ))
    return definition
end

function _parse_file_input(entry::AbstractDict, bindings::NamedTuple, index::Int)
    context = "inputs[$index]"
    fields = ("name", "destination", "binding")
    _require_table_keys(entry, fields, fields, context)
    name = _file_identifier(entry["name"], "$context.name")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    values = _require_binding(bindings, entry["binding"], "$context.binding")
    return graph_input(name, destination, values)
end

function _parse_file_output(entry::AbstractDict, index::Int)
    context = "outputs[$index]"
    fields = ("name", "source")
    _require_table_keys(entry, fields, fields, context)
    name = _file_identifier(entry["name"], "$context.name")
    source = _file_endpoint(entry["source"], "$context.source")
    return graph_output(name, source)
end

function _parse_file_link(entry::AbstractDict, index::Int)
    context = "links[$index]"
    fields = ("source", "destination")
    _require_table_keys(entry, fields, fields, context)
    source = _file_endpoint(entry["source"], "$context.source")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    return link(source, destination)
end

function _parse_file_delayed_link(
    entry::AbstractDict,
    bindings::NamedTuple,
    index::Int,
)
    context = "delayed_links[$index]"
    fields = ("source", "destination", "initial")
    _require_table_keys(entry, fields, fields, context)
    source = _file_endpoint(entry["source"], "$context.source")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    initial = _require_binding(bindings, entry["initial"], "$context.initial")
    return delayed_link(source, destination, initial)
end

function _parse_file_parameter(
    entry::AbstractDict,
    bindings::NamedTuple,
    index::Int,
)
    context = "parameters[$index]"
    fields = ("destination", "binding")
    _require_table_keys(entry, fields, fields, context)
    destination = _file_endpoint(entry["destination"], "$context.destination")
    values = _require_binding(bindings, entry["binding"], "$context.binding")
    return sparse_parameter(destination, values)
end

function _parse_entries(entries, parse_entry)
    return Tuple(
        parse_entry(entry, Int(index))
        for (index, entry) in pairs(entries)
    )
end

function _read_algorithm_graph(path)
    throw(AlgorithmGraphError(
        "YAML graph file loading requires YAML.jl; load it with `using YAML`",
    ))
end

"""
    load_algorithm_graph(path; node_types=builtin_graph_node_types(),
                         bindings=NamedTuple())
    load_algorithm_graph(root::AbstractDict;
                         node_types=builtin_graph_node_types(),
                         bindings=NamedTuple())

Load a versioned, declarative graph and compile it into the same concrete
[`AlgorithmGraphDefinition`](@ref) used by direct Julia construction. The path
method parses YAML after `using YAML` activates the optional file extension;
the dictionary method accepts an already-parsed root, such
as one decoded by an application-owned interchange layer.
`node_types` is an explicit immutable named tuple of cold factories. Detector
factories accept `(name, config, props; detector)`. Other factories accept
`(name, config, props)` or `(name, config, props, bindings)`; the latter is
selected when applicable and can resolve explicit calibration binding names.
[`graph_node_bindings`](@ref) declares the construction calibrations required
by such a factory without loading them.
`bindings` supplies caller-owned frame inputs, delayed-link initial values,
sparse parameters, and cold node calibrations by name; the graph dictionary
never embeds or implicitly loads calibration arrays.
"""
function load_algorithm_graph(
    path::AbstractString;
    node_types::NamedTuple=builtin_graph_node_types(),
    bindings::NamedTuple=NamedTuple(),
)
    root = try
        _read_algorithm_graph(path)
    catch error
        throw(AlgorithmGraphError(
            "could not parse algorithm graph '$path': $(sprint(showerror, error))",
        ))
    end
    return _load_algorithm_graph(
        root;
        node_types,
        bindings,
        context="algorithm graph '$path'",
    )
end

function load_algorithm_graph(
    root::AbstractDict;
    node_types::NamedTuple=builtin_graph_node_types(),
    bindings::NamedTuple=NamedTuple(),
)
    return _load_algorithm_graph(
        root;
        node_types,
        bindings,
        context="algorithm graph dictionary",
    )
end

function _load_algorithm_graph(
    root::AbstractDict;
    node_types::NamedTuple,
    bindings::NamedTuple,
    context::AbstractString,
)
    _require_table_keys(
        root,
        (
            "schema_version",
            "name",
            "nodes",
            "inputs",
            "outputs",
            "links",
            "delayed_links",
            "parameters",
        ),
        ("schema_version", "name", "nodes"),
        context,
    )
    version = _file_integer(root["schema_version"], "schema_version")
    version == ALGORITHM_GRAPH_SCHEMA_VERSION || throw(AlgorithmGraphError(
        "unsupported algorithm graph schema_version $version; expected " *
        "$ALGORITHM_GRAPH_SCHEMA_VERSION",
    ))
    name = _file_identifier(root["name"], "name")

    node_entries = _table_array_or_empty(root, "nodes", context)
    input_entries = _table_array_or_empty(root, "inputs", context)
    output_entries = _table_array_or_empty(root, "outputs", context)
    link_entries = _table_array_or_empty(root, "links", context)
    delayed_entries = _table_array_or_empty(root, "delayed_links", context)
    parameter_entries = _table_array_or_empty(root, "parameters", context)

    nodes = _parse_entries(
        node_entries,
        (entry, index) -> _parse_file_node(entry, node_types, bindings, index),
    )
    inputs = _parse_entries(
        input_entries,
        (entry, index) -> _parse_file_input(entry, bindings, index),
    )
    outputs = _parse_entries(output_entries, _parse_file_output)
    links = _parse_entries(link_entries, _parse_file_link)
    delayed_links = _parse_entries(
        delayed_entries,
        (entry, index) -> _parse_file_delayed_link(entry, bindings, index),
    )
    parameters = _parse_entries(
        parameter_entries,
        (entry, index) -> _parse_file_parameter(entry, bindings, index),
    )
    return algorithm_graph(
        nodes;
        name=name,
        inputs=inputs,
        outputs=outputs,
        links=links,
        delayed_links=delayed_links,
        parameters=parameters,
    )
end

function _load_algorithm_graph(
    root;
    node_types::NamedTuple,
    bindings::NamedTuple,
    context::AbstractString,
)
    throw(_file_error(context, "expected a root mapping, got $(typeof(root))"))
end

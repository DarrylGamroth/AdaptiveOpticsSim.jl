const ALGORITHM_GRAPH_TOML_SCHEMA_VERSION = 1

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

_normalize_toml_value(value::Union{Bool,Integer,AbstractFloat,AbstractString}) =
    value

function _normalize_toml_value(values::AbstractVector)
    builder = Vector{Any}(undef, length(values))
    for index in eachindex(values)
        builder[index] = _normalize_toml_value(values[index])
    end
    return Tuple(builder)
end


function _normalize_toml_value(table::AbstractDict)
    names = sort!(collect(String, keys(table)))
    symbols = Tuple(Symbol(name) for name in names)
    values = Tuple(_normalize_toml_value(table[name]) for name in names)
    return NamedTuple{symbols}(values)
end

function _normalize_toml_value(value)
    throw(AlgorithmGraphError(
        "unsupported TOML value type $(typeof(value)) in graph config or props",
    ))
end

function _table_or_empty(table::AbstractDict, key::String, context::AbstractString)
    value = get(table, key, Dict{String,Any}())
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
    value = get(table, key, Any[])
    value isa AbstractVector || throw(_file_error(
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

function _discrete_integrator_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :extent,
        :input_schema,
        :output_schema,
        :sample_period_s,
    )
    prop_fields = (:gain, :tau_s)
    _require_named_fields(
        config,
        config_fields,
        config_fields,
        "node '$name' config",
    )
    _require_named_fields(
        props,
        prop_fields,
        prop_fields,
        "node '$name' props",
    )
    return discrete_integrator_node(
        name;
        extent=_file_integer(config.extent, "node '$name' config.extent"),
        sample_period_s=_file_real(
            config.sample_period_s,
            "node '$name' config.sample_period_s",
        ),
        input_schema=_file_string(
            config.input_schema,
            "node '$name' config.input_schema",
        ),
        output_schema=_file_string(
            config.output_schema,
            "node '$name' config.output_schema",
        ),
        gain=_file_real(props.gain, "node '$name' props.gain"),
        tau_s=_file_real(props.tau_s, "node '$name' props.tau_s"),
        T=Float32,
    )
end

function _closed_loop_correction_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :constraint_feedback_schema,
        :controller_state_schema,
        :correction_schema,
        :extent,
        :residual_error_schema,
    )
    prop_fields = (:anti_windup_gain, :gain, :pole)
    context = "node '$name' config"
    props_context = "node '$name' props"
    _require_named_fields(config, config_fields, config_fields, context)
    _require_named_fields(props, prop_fields, prop_fields, props_context)
    return closed_loop_correction_node(
        name;
        extent=_file_integer(config.extent, "$context.extent"),
        residual_error_schema=_file_string(
            config.residual_error_schema,
            "$context.residual_error_schema",
        ),
        constraint_feedback_schema=_file_string(
            config.constraint_feedback_schema,
            "$context.constraint_feedback_schema",
        ),
        correction_schema=_file_string(
            config.correction_schema,
            "$context.correction_schema",
        ),
        controller_state_schema=_file_string(
            config.controller_state_schema,
            "$context.controller_state_schema",
        ),
        gain=_file_real(props.gain, "$props_context.gain"),
        pole=_file_real(props.pole, "$props_context.pole"),
        anti_windup_gain=_file_real(
            props.anti_windup_gain,
            "$props_context.anti_windup_gain",
        ),
        T=Float32,
    )
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
        :light_ratio,
        :modulation,
        :modulation_points,
        :n_pix_edge,
        :n_pix_separation,
        :opd_schema,
        :photon_rate_schema,
        :psf_centering,
        :pupil_reflectivity,
        :pupil_samples,
        :resolution,
        :source_band,
        :source_magnitude,
        :source_photon_irradiance_m2_s,
        :source_position_angle_deg,
        :source_separation_arcsec,
        :source_wavelength_m,
        :telescope_diameter_m,
        :threshold,
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
    threshold = hasproperty(config, :threshold) ?
        _file_real(config.threshold, "$context.threshold") : 0.1
    modulation = hasproperty(config, :modulation) ?
        _file_real(config.modulation, "$context.modulation") : 2.0
    modulation_points = hasproperty(config, :modulation_points) ?
        _file_integer(
            config.modulation_points,
            "$context.modulation_points",
        ) : nothing
    light_ratio = hasproperty(config, :light_ratio) ?
        _file_real(config.light_ratio, "$context.light_ratio") : 0.0
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
        threshold,
        modulation,
        modulation_points,
        light_ratio,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge,
        binning,
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

function _shack_hartmann_centroid_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :calibration_signature,
        :calibration_wavelength_m,
        :centroid_response,
        :centroid_cutoff_fraction,
        :frame_schema,
        :n_lenslets,
        :n_pix_subap,
        :reference_signal_schema,
        :resolution,
        :slopes_schema,
        :telescope_diameter_m,
        :valid_subapertures_schema,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, config_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    return shack_hartmann_centroid_node(
        name;
        resolution=_file_integer(config.resolution, "$context.resolution"),
        telescope_diameter_m=_file_real(
            config.telescope_diameter_m,
            "$context.telescope_diameter_m",
        ),
        n_lenslets=_file_integer(
            config.n_lenslets,
            "$context.n_lenslets",
        ),
        n_pix_subap=_file_integer(
            config.n_pix_subap,
            "$context.n_pix_subap",
        ),
        centroid_cutoff_fraction=_file_real(
            config.centroid_cutoff_fraction,
            "$context.centroid_cutoff_fraction",
        ),
        centroid_response=_file_real(
            config.centroid_response,
            "$context.centroid_response",
        ),
        calibration_wavelength_m=_file_real(
            config.calibration_wavelength_m,
            "$context.calibration_wavelength_m",
        ),
        calibration_signature=_file_integer(
            config.calibration_signature,
            "$context.calibration_signature",
        ),
        frame_schema=_file_string(
            config.frame_schema,
            "$context.frame_schema",
        ),
        slopes_schema=_file_string(
            config.slopes_schema,
            "$context.slopes_schema",
        ),
        valid_subapertures_schema=_file_string(
            config.valid_subapertures_schema,
            "$context.valid_subapertures_schema",
        ),
        reference_signal_schema=_file_string(
            config.reference_signal_schema,
            "$context.reference_signal_schema",
        ),
        T=Float32,
    )
end

function _shack_hartmann_slope_selection_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :full_slopes_schema,
        :lenslet_order_schema,
        :n_lenslets,
        :selected_lenslet_count,
        :selected_slopes_schema,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, config_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    return shack_hartmann_slope_selection_node(
        name;
        n_lenslets=_file_integer(
            config.n_lenslets,
            "$context.n_lenslets",
        ),
        selected_lenslet_count=_file_integer(
            config.selected_lenslet_count,
            "$context.selected_lenslet_count",
        ),
        full_slopes_schema=_file_string(
            config.full_slopes_schema,
            "$context.full_slopes_schema",
        ),
        selected_slopes_schema=_file_string(
            config.selected_slopes_schema,
            "$context.selected_slopes_schema",
        ),
        lenslet_order_schema=_file_string(
            config.lenslet_order_schema,
            "$context.lenslet_order_schema",
        ),
        T=Float32,
    )
end

function _control_matrix_reconstruction_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :control_matrix_schema,
        :reconstructed_count,
        :reconstructed_schema,
        :slope_count,
        :slopes_schema,
    )
    context = "node '$name' config"
    _require_named_fields(config, config_fields, config_fields, context)
    _require_named_fields(props, (), (), "node '$name' props")
    return control_matrix_reconstruction_node(
        name;
        slope_count=_file_integer(
            config.slope_count,
            "$context.slope_count",
        ),
        reconstructed_count=_file_integer(
            config.reconstructed_count,
            "$context.reconstructed_count",
        ),
        slopes_schema=_file_string(
            config.slopes_schema,
            "$context.slopes_schema",
        ),
        reconstructed_schema=_file_string(
            config.reconstructed_schema,
            "$context.reconstructed_schema",
        ),
        control_matrix_schema=_file_string(
            config.control_matrix_schema,
            "$context.control_matrix_schema",
        ),
        T=Float32,
    )
end

"""
    builtin_graph_node_types()

Return the immutable node-type map accepted by the maintained TOML loader.
Optional packages may merge additional cold factories into this named tuple
and pass the result explicitly to [`load_algorithm_graph`](@ref). No mutable
global registry or dynamic Julia evaluation is used.
"""
function builtin_graph_node_types()
    return (
        ccd_detector_acquisition_f32=_ccd_detector_acquisition_f32_node,
        closed_loop_correction_f32=_closed_loop_correction_f32_node,
        control_matrix_reconstruction_f32=
            _control_matrix_reconstruction_f32_node,
        discrete_integrator_f32=_discrete_integrator_f32_node,
        emccd_detector_acquisition_f32=
            _emccd_detector_acquisition_f32_node,
        modal_opd_expansion_f32=_modal_opd_expansion_f32_node,
        pyramid_rate_f32=_pyramid_rate_f32_node,
        shack_hartmann_centroid_f32=_shack_hartmann_centroid_f32_node,
        shack_hartmann_rate_f32=_shack_hartmann_rate_f32_node,
        shack_hartmann_slope_selection_f32=
            _shack_hartmann_slope_selection_f32_node,
    )
end

function _parse_file_node(entry::AbstractDict, node_types::NamedTuple, index::Int)
    context = "nodes[$index]"
    _require_table_keys(
        entry,
        ("name", "type", "config", "props"),
        ("name", "type", "config"),
        context,
    )
    name = _file_identifier(entry["name"], "$context.name")
    node_type = _file_identifier(entry["type"], "$context.type")
    hasproperty(node_types, node_type) || throw(_file_error(
        "$context.type",
        "unknown graph node type '$node_type'",
    ))
    config = _normalize_toml_value(
        _table_or_empty(entry, "config", context),
    )
    props = _normalize_toml_value(
        _table_or_empty(entry, "props", context),
    )
    factory = getproperty(node_types, node_type)
    definition = factory(name, config, props)
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

function _parse_entries(entries::AbstractVector, parse_entry)
    builder = Vector{Any}(undef, length(entries))
    for (index, entry) in pairs(entries)
        builder[index] = parse_entry(entry, Int(index))
    end
    return Tuple(builder)
end

"""
    load_algorithm_graph(path; node_types=builtin_graph_node_types(),
                         bindings=NamedTuple())

Load a versioned, declarative TOML graph and compile it into the same concrete
[`AlgorithmGraphDefinition`](@ref) used by direct Julia construction.
`node_types` is an explicit immutable named tuple of cold factories. `bindings`
supplies caller-owned frame inputs, delayed-link initial values, and large
sparse parameters by name; the TOML file never embeds or implicitly loads
calibration arrays.
"""
function load_algorithm_graph(
    path::AbstractString;
    node_types::NamedTuple=builtin_graph_node_types(),
    bindings::NamedTuple=NamedTuple(),
)
    root = try
        TOML.parsefile(path)
    catch error
        throw(AlgorithmGraphError(
            "could not parse algorithm graph '$path': $(sprint(showerror, error))",
        ))
    end
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
        "algorithm graph '$path'",
    )
    version = _file_integer(root["schema_version"], "schema_version")
    version == ALGORITHM_GRAPH_TOML_SCHEMA_VERSION || throw(AlgorithmGraphError(
        "unsupported algorithm graph schema_version $version; expected " *
        "$ALGORITHM_GRAPH_TOML_SCHEMA_VERSION",
    ))
    name = _file_identifier(root["name"], "name")

    node_entries = _table_array_or_empty(root, "nodes", "algorithm graph '$path'")
    input_entries = _table_array_or_empty(root, "inputs", "algorithm graph '$path'")
    output_entries = _table_array_or_empty(root, "outputs", "algorithm graph '$path'")
    link_entries = _table_array_or_empty(root, "links", "algorithm graph '$path'")
    delayed_entries = _table_array_or_empty(
        root,
        "delayed_links",
        "algorithm graph '$path'",
    )
    parameter_entries = _table_array_or_empty(
        root,
        "parameters",
        "algorithm graph '$path'",
    )

    nodes = _parse_entries(
        node_entries,
        (entry, index) -> _parse_file_node(entry, node_types, index),
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

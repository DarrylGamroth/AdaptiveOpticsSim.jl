struct DiscreteIntegratorNode{T<:AbstractFloat} end

"""Construction and graph-rebuild values for one discrete-integrator node."""
struct DiscreteIntegratorNodeConfig{T<:AbstractFloat}
    extent::Int
    sample_period_s::T
    input_schema::String
    output_schema::String

    function DiscreteIntegratorNodeConfig{T}(
        extent::Integer,
        sample_period_s::Real,
        input_schema::AbstractString,
        output_schema::AbstractString,
    ) where {T<:AbstractFloat}
        extent > 0 || throw(AlgorithmGraphError(
            "discrete-integrator extent must be greater than zero",
        ))
        period = T(sample_period_s)
        isfinite(period) && period > zero(T) || throw(AlgorithmGraphError(
            "discrete-integrator sample_period_s must be finite and greater than zero",
        ))
        isempty(input_schema) && throw(AlgorithmGraphError(
            "discrete-integrator input_schema must not be empty",
        ))
        isempty(output_schema) && throw(AlgorithmGraphError(
            "discrete-integrator output_schema must not be empty",
        ))
        return new{T}(
            Int(extent),
            period,
            String(input_schema),
            String(output_schema),
        )
    end
end

"""Initial scalar props for one discrete-integrator node."""
struct DiscreteIntegratorNodeProps{T<:AbstractFloat}
    gain::T
    tau_s::T

    function DiscreteIntegratorNodeProps{T}(
        gain::Real,
        tau_s::Real,
    ) where {T<:AbstractFloat}
        converted_gain = T(gain)
        converted_tau_s = T(tau_s)
        try
            DiscreteIntegratorPlan(converted_gain, converted_tau_s)
        catch error
            throw(AlgorithmGraphError(sprint(showerror, error)))
        end
        return new{T}(converted_gain, converted_tau_s)
    end
end

"""
    discrete_integrator_node(name; extent, sample_period_s, input_schema,
                             output_schema, gain=0.3, tau_s=0.02,
                             T=Float32)

Declare one fixed-step graph node around the native AOS discrete-integrator
implementation. The node owns its plan, persistent controller state,
replaceable workspace, and exact graph-buffer bindings after preparation.
"""
function discrete_integrator_node(
    name::Symbol;
    extent::Integer,
    sample_period_s::Real,
    input_schema::AbstractString,
    output_schema::AbstractString,
    gain::Real=0.3,
    tau_s::Real=0.02,
    T::Type{<:AbstractFloat}=Float32,
)
    config = DiscreteIntegratorNodeConfig{T}(
        extent,
        sample_period_s,
        input_schema,
        output_schema,
    )
    props = DiscreteIntegratorNodeProps{T}(gain, tau_s)
    return algorithm_node(name, DiscreteIntegratorNode{T}, config; props=props)
end

function graph_node_ports(
    ::Type{DiscreteIntegratorNode{T}},
    config::DiscreteIntegratorNodeConfig{T},
) where {T}
    format(name, direction, schema) = graph_port_contract(
        name,
        direction,
        :data,
        T,
        (config.extent,),
        schema,
        :column_major,
    )
    return (
        format(:input, :input, config.input_schema),
        format(:output, :output, config.output_schema),
    )
end

struct _DiscreteIntegratorStepPlan{
    T<:AbstractFloat,
    P<:DiscreteIntegratorPlan{T},
}
    controller::P
    sample_period_s::T
end

struct _DiscreteIntegratorOwner{Plan,State,Workspace,Input,Output}
    plan::Plan
    state::State
    workspace::Workspace
    input::Input
    output::Output
end

function prepare_graph_node(
    ::Type{DiscreteIntegratorNode{T}},
    config::DiscreteIntegratorNodeConfig{T},
    props::DiscreteIntegratorNodeProps{T},
    inputs::NamedTuple{(:input,)},
    outputs::NamedTuple{(:output,)},
    ::NamedTuple{()},
    target,
) where {T}
    plan = _DiscreteIntegratorStepPlan(
        DiscreteIntegratorPlan(props.gain, props.tau_s),
        config.sample_period_s,
    )
    state = DiscreteIntegratorState(outputs.output)
    workspace = DiscreteIntegratorWorkspace(outputs.output)
    return _DiscreteIntegratorOwner(
        plan,
        state,
        workspace,
        inputs.input,
        outputs.output,
    )
end

@inline function step_graph_node!(owner::_DiscreteIntegratorOwner)
    update!(
        owner.output,
        owner.state,
        owner.workspace,
        owner.plan.controller,
        owner.input,
        owner.plan.sample_period_s,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_DiscreteIntegratorOwner)
    reset_controller!(owner.state, owner.workspace)
    return nothing
end

struct ModalOPDExpansionNode{T<:AbstractFloat} end

"""Construction and graph-rebuild values for one modal OPD expansion node."""
struct ModalOPDExpansionNodeConfig
    pupil_rows::Int
    pupil_columns::Int
    mode_count::Int
    coefficients_schema::String
    opd_schema::String
    basis_schema::String
    pupil_support_schema::String

    function ModalOPDExpansionNodeConfig(
        pupil_rows::Integer,
        pupil_columns::Integer,
        mode_count::Integer,
        coefficients_schema::AbstractString,
        opd_schema::AbstractString,
        basis_schema::AbstractString,
        pupil_support_schema::AbstractString,
    )
        pupil_rows > 0 || throw(AlgorithmGraphError(
            "modal OPD pupil_rows must be greater than zero",
        ))
        pupil_columns > 0 || throw(AlgorithmGraphError(
            "modal OPD pupil_columns must be greater than zero",
        ))
        mode_count > 0 || throw(AlgorithmGraphError(
            "modal OPD mode_count must be greater than zero",
        ))
        schemas = (
            coefficients_schema,
            opd_schema,
            basis_schema,
            pupil_support_schema,
        )
        all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
            "modal OPD scientific schemas must not be empty",
        ))
        return new(
            Int(pupil_rows),
            Int(pupil_columns),
            Int(mode_count),
            String(coefficients_schema),
            String(opd_schema),
            String(basis_schema),
            String(pupil_support_schema),
        )
    end
end

"""
    modal_opd_expansion_node(name; pupil_rows, pupil_columns, mode_count,
                             coefficients_schema, opd_schema, basis_schema,
                             pupil_support_schema, T=Float32)

Declare one stateless graph node around the native AOS modal OPD expansion.
The graph supplies the modal basis and pupil support as explicit startup sparse
parameters.
"""
function modal_opd_expansion_node(
    name::Symbol;
    pupil_rows::Integer,
    pupil_columns::Integer,
    mode_count::Integer,
    coefficients_schema::AbstractString,
    opd_schema::AbstractString,
    basis_schema::AbstractString,
    pupil_support_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ModalOPDExpansionNodeConfig(
        pupil_rows,
        pupil_columns,
        mode_count,
        coefficients_schema,
        opd_schema,
        basis_schema,
        pupil_support_schema,
    )
    return algorithm_node(
        name,
        ModalOPDExpansionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ModalOPDExpansionNode{T}},
    config::ModalOPDExpansionNodeConfig,
) where {T}
    return (
        graph_port_contract(
            :coefficients,
            :input,
            :data,
            T,
            (config.mode_count,),
            config.coefficients_schema,
            :column_major,
        ),
        graph_port_contract(
            :opd,
            :output,
            :data,
            T,
            (config.pupil_rows, config.pupil_columns),
            config.opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :basis,
            :input,
            :parameter,
            T,
            (config.pupil_rows, config.pupil_columns, config.mode_count),
            config.basis_schema,
            :column_major,
        ),
        graph_port_contract(
            :pupil_support,
            :input,
            :parameter,
            Bool,
            (config.pupil_rows, config.pupil_columns),
            config.pupil_support_schema,
            :column_major,
        ),
    )
end

struct _ModalOPDExpansionOwner{Plan,Coefficients,OPD}
    plan::Plan
    coefficients::Coefficients
    opd::OPD
end

function prepare_graph_node(
    ::Type{ModalOPDExpansionNode{T}},
    ::ModalOPDExpansionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:coefficients,)},
    outputs::NamedTuple{(:opd,)},
    parameters::NamedTuple{(:basis,:pupil_support)},
    target,
) where {T}
    plan = ModalOPDExpansionPlan(parameters.basis, parameters.pupil_support)
    return _ModalOPDExpansionOwner(plan, inputs.coefficients, outputs.opd)
end

@inline function step_graph_node!(owner::_ModalOPDExpansionOwner)
    combine_basis!(owner.opd, owner.plan, owner.coefficients)
    return nothing
end

@inline reset_graph_node!(::_ModalOPDExpansionOwner) = nothing

struct ShackHartmannRateNode{T<:AbstractFloat} end

"""Construction values for one diffractive Shack–Hartmann photon-rate node."""
struct ShackHartmannRateNodeConfig{T<:AbstractFloat,TD,S}
    telescope::TD
    source::S
    n_lenslets::Int
    n_pix_subap::Int
    diffraction_padding::Int
    pixel_scale_arcsec::T
    valid_subaperture_threshold::T
    threshold_convolution::T
    half_pixel_shift::Bool
    shannon_sampling::Bool
    opd_schema::String
    photon_rate_schema::String
end

function _shack_hartmann_rate_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    diffraction_padding::Integer,
    pixel_scale_arcsec::Real,
    valid_subaperture_threshold::Real,
    threshold_convolution::Real,
    half_pixel_shift::Bool,
    shannon_sampling::Bool,
    source_band::Symbol,
    source_magnitude::Real,
    source_wavelength_m::Real,
    source_photon_irradiance_m2_s::Union{Nothing,Real},
    source_separation_arcsec::Real,
    source_position_angle_deg::Real,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
) where {T<:AbstractFloat}
    lenslets = Int(n_lenslets)
    pixels = Int(n_pix_subap)
    padding = Int(diffraction_padding)
    pupil_resolution = Int(resolution)
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann n_lenslets must be greater than zero",
    ))
    pupil_resolution > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann resolution must be greater than zero",
    ))
    pupil_resolution % lenslets == 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann resolution must be divisible by n_lenslets",
    ))
    pixels > 0 && iseven(pixels) || throw(AlgorithmGraphError(
        "Shack–Hartmann n_pix_subap must be positive and even",
    ))
    padding > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann diffraction_padding must be greater than zero",
    ))
    isempty(opd_schema) && throw(AlgorithmGraphError(
        "Shack–Hartmann opd_schema must not be empty",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "Shack–Hartmann photon_rate_schema must not be empty",
    ))

    pixel_scale = T(pixel_scale_arcsec)
    valid_threshold = T(valid_subaperture_threshold)
    convolution_threshold = T(threshold_convolution)
    isfinite(pixel_scale) && pixel_scale > zero(T) || throw(
        AlgorithmGraphError(
            "Shack–Hartmann pixel_scale_arcsec must be finite and positive",
        ))
    isfinite(valid_threshold) && zero(T) <= valid_threshold <= one(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann valid_subaperture_threshold must lie in [0, 1]",
        ))
    isfinite(convolution_threshold) &&
        zero(T) <= convolution_threshold <= one(T) || throw(
            AlgorithmGraphError(
                "Shack–Hartmann threshold_convolution must lie in [0, 1]",
            ))

    telescope, source = try
        definition = TelescopeDefinition(
            resolution=pupil_resolution,
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=aperture_revision,
            T=T,
        )
        guide_source = Source(
            band=source_band,
            magnitude=source_magnitude,
            coordinates=(source_separation_arcsec, source_position_angle_deg),
            wavelength=source_wavelength_m,
            photon_irradiance=source_photon_irradiance_m2_s,
            T=T,
        )
        photon_irradiance(guide_source)
        (definition, guide_source)
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end

    return ShackHartmannRateNodeConfig(
        telescope,
        source,
        lenslets,
        pixels,
        padding,
        pixel_scale,
        valid_threshold,
        convolution_threshold,
        half_pixel_shift,
        shannon_sampling,
        String(opd_schema),
        String(photon_rate_schema),
    )
end

"""
    shack_hartmann_rate_node(name; ...)

Declare one complete-frame diffractive Shack–Hartmann optics node. The node
consumes a pupil-plane OPD map and writes a cell-integrated detector-plane
photon-rate mosaic. Detector exposure/noise and slope estimation remain
separate graph operations.
"""
function shack_hartmann_rate_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    pixel_scale_arcsec::Real,
    source_wavelength_m::Real,
    source_photon_irradiance_m2_s::Union{Nothing,Real}=nothing,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    diffraction_padding::Integer=2,
    valid_subaperture_threshold::Real=0.1,
    threshold_convolution::Real=0.05,
    half_pixel_shift::Bool=false,
    shannon_sampling::Bool=true,
    source_band::Symbol=:custom,
    source_magnitude::Real=0,
    source_separation_arcsec::Real=0,
    source_position_angle_deg::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _shack_hartmann_rate_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        n_lenslets,
        n_pix_subap,
        diffraction_padding,
        pixel_scale_arcsec,
        valid_subaperture_threshold,
        threshold_convolution,
        half_pixel_shift,
        shannon_sampling,
        source_band,
        source_magnitude,
        source_wavelength_m,
        source_photon_irradiance_m2_s,
        source_separation_arcsec,
        source_position_angle_deg,
        opd_schema,
        photon_rate_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannRateNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannRateNode{T}},
    config::ShackHartmannRateNodeConfig{T},
) where {T}
    output_extent = config.n_lenslets * config.n_pix_subap
    return (
        graph_port_contract(
            :opd,
            :input,
            :data,
            T,
            (config.telescope.resolution, config.telescope.resolution),
            config.opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :photon_rate,
            :output,
            :data,
            T,
            (output_extent, output_extent),
            config.photon_rate_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannRateOwner{P}
    prepared::P
end

function prepare_graph_node(
    ::Type{ShackHartmannRateNode{T}},
    config::ShackHartmannRateNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:opd,)},
    outputs::NamedTuple{(:photon_rate,)},
    ::NamedTuple{()},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    pupil = PupilFunction(telescope, inputs.opd)
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=config.n_lenslets,
        threshold=config.valid_subaperture_threshold,
        threshold_convolution=config.threshold_convolution,
        half_pixel_shift=config.half_pixel_shift,
        diffraction_padding=config.diffraction_padding,
        pixel_scale_arcsec=config.pixel_scale_arcsec,
        n_pix_subap=config.n_pix_subap,
        shannon_sampling=config.shannon_sampling,
        mode=Diffractive(),
        T=T,
        backend=compute_device_backend(target),
    )
    optics = shack_hartmann_optics(sensor, config.source)
    photon_rate = shack_hartmann_rate_map(
        optics,
        pupil,
        outputs.photon_rate,
    )
    prepared = prepare_wfs_optics(optics, pupil, photon_rate)
    return _ShackHartmannRateOwner(prepared)
end

@inline function step_graph_node!(owner::_ShackHartmannRateOwner)
    prepared = owner.prepared
    form_wfs_optical_products!(
        prepared.output,
        prepared.input,
        prepared,
    )
    return nothing
end

@inline reset_graph_node!(::_ShackHartmannRateOwner) = nothing

struct ShackHartmannCentroidNode{T<:AbstractFloat} end

"""Construction values for one calibrated Shack–Hartmann centroid node."""
struct ShackHartmannCentroidNodeConfig{T<:AbstractFloat,TD}
    telescope::TD
    n_lenslets::Int
    n_pix_subap::Int
    centroid_cutoff_fraction::T
    centroid_response::T
    calibration_wavelength_m::T
    calibration_signature::UInt
    frame_schema::String
    slopes_schema::String
    valid_subapertures_schema::String
    reference_signal_schema::String
end

function _shack_hartmann_centroid_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    centroid_cutoff_fraction::Real,
    centroid_response::Real,
    calibration_wavelength_m::Real,
    calibration_signature::Integer,
    frame_schema::AbstractString,
    slopes_schema::AbstractString,
    valid_subapertures_schema::AbstractString,
    reference_signal_schema::AbstractString,
) where {T<:AbstractFloat}
    pupil_resolution = Int(resolution)
    lenslets = Int(n_lenslets)
    pixels = Int(n_pix_subap)
    pupil_resolution > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid resolution must be greater than zero",
    ))
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid n_lenslets must be greater than zero",
    ))
    pupil_resolution % lenslets == 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid resolution must be divisible by n_lenslets",
    ))
    pixels > 0 && iseven(pixels) || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid n_pix_subap must be positive and even",
    ))

    cutoff_fraction = T(centroid_cutoff_fraction)
    response = T(centroid_response)
    wavelength = T(calibration_wavelength_m)
    isfinite(cutoff_fraction) &&
        zero(T) <= cutoff_fraction <= one(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann centroid_cutoff_fraction must lie in [0, 1]",
        ))
    isfinite(response) && response != zero(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann centroid_response must be finite and nonzero",
        ))
    isfinite(wavelength) && wavelength > zero(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann calibration_wavelength_m must be finite and positive",
        ))
    calibration_signature >= 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann calibration_signature must be nonnegative",
    ))
    calibration_signature <= typemax(UInt) || throw(AlgorithmGraphError(
        "Shack–Hartmann calibration_signature exceeds UInt",
    ))
    schemas = (
        frame_schema,
        slopes_schema,
        valid_subapertures_schema,
        reference_signal_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid scientific schemas must not be empty",
    ))

    telescope = try
        TelescopeDefinition(
            resolution=pupil_resolution,
            diameter=telescope_diameter_m,
            revision=0,
            T=T,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return ShackHartmannCentroidNodeConfig(
        telescope,
        lenslets,
        pixels,
        cutoff_fraction,
        response,
        wavelength,
        UInt(calibration_signature),
        String(frame_schema),
        String(slopes_schema),
        String(valid_subapertures_schema),
        String(reference_signal_schema),
    )
end

"""
    shack_hartmann_centroid_node(name; ...)

Declare one calibrated, complete-frame Shack–Hartmann centroid operation. The
node consumes a lenslet-mosaic detector frame and writes the full canonical
`[axis 1; axis 2]` centroid-slope vector. Its valid-subaperture mask and
reference signal are explicit startup sparse parameters; compact slope
selection remains a separate graph operation.
"""
function shack_hartmann_centroid_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    centroid_cutoff_fraction::Real,
    centroid_response::Real,
    calibration_wavelength_m::Real,
    calibration_signature::Integer,
    frame_schema::AbstractString,
    slopes_schema::AbstractString,
    valid_subapertures_schema::AbstractString,
    reference_signal_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _shack_hartmann_centroid_config(
        T;
        resolution,
        telescope_diameter_m,
        n_lenslets,
        n_pix_subap,
        centroid_cutoff_fraction,
        centroid_response,
        calibration_wavelength_m,
        calibration_signature,
        frame_schema,
        slopes_schema,
        valid_subapertures_schema,
        reference_signal_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannCentroidNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannCentroidNode{T}},
    config::ShackHartmannCentroidNodeConfig{T},
) where {T}
    lenslet_count = config.n_lenslets * config.n_lenslets
    frame_extent = config.n_lenslets * config.n_pix_subap
    return (
        graph_port_contract(
            :frame,
            :input,
            :data,
            T,
            (frame_extent, frame_extent),
            config.frame_schema,
            :column_major,
        ),
        graph_port_contract(
            :slopes,
            :output,
            :data,
            T,
            (2 * lenslet_count,),
            config.slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :valid_subapertures,
            :input,
            :parameter,
            Bool,
            (config.n_lenslets, config.n_lenslets),
            config.valid_subapertures_schema,
            :column_major,
        ),
        graph_port_contract(
            :reference_signal,
            :input,
            :parameter,
            T,
            (lenslet_count, 2),
            config.reference_signal_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannCentroidOwner{P,O,M}
    prepared::P
    observation::O
    measurement::M
end

function prepare_graph_node(
    ::Type{ShackHartmannCentroidNode{T}},
    config::ShackHartmannCentroidNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:frame,)},
    outputs::NamedTuple{(:slopes,)},
    parameters::NamedTuple{(:valid_subapertures,:reference_signal)},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=config.n_lenslets,
        threshold_cog=config.centroid_cutoff_fraction,
        diffraction_padding=1,
        n_pix_subap=config.n_pix_subap,
        mode=Diffractive(),
        T=T,
        backend=compute_device_backend(target),
    )
    set_valid_subapertures!(sensor, parameters.valid_subapertures)
    set_subaperture_calibration!(
        subaperture_calibration(sensor),
        parameters.reference_signal;
        centroid_response=config.centroid_response,
        output_units=:pixel,
        wavelength=config.calibration_wavelength_m,
        signature=config.calibration_signature,
    )
    observation = WFSObservation(
        inputs.frame;
        units=:electron_count,
        layout=:lenslet_mosaic,
    )
    measurement = WFSMeasurement(
        outputs.slopes;
        units=:pixel,
        kind=:centroid_slopes,
    )
    prepared = prepare_wfs_estimation(sensor, observation, measurement)
    return _ShackHartmannCentroidOwner(
        prepared,
        observation,
        measurement,
    )
end

@inline function step_graph_node!(owner::_ShackHartmannCentroidOwner)
    estimate_wfs_measurement!(
        owner.measurement,
        owner.observation,
        owner.prepared,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_ShackHartmannCentroidOwner)
    fill!(owner.measurement.storage, zero(eltype(owner.measurement.storage)))
    return nothing
end

struct ShackHartmannSlopeSelectionNode{T<:AbstractFloat} end

"""Construction values for one ordered Shack–Hartmann slope selector."""
struct ShackHartmannSlopeSelectionNodeConfig
    n_lenslets::Int
    selected_lenslet_count::Int
    full_slopes_schema::String
    selected_slopes_schema::String
    lenslet_order_schema::String
end

function ShackHartmannSlopeSelectionNodeConfig(
    n_lenslets::Integer,
    selected_lenslet_count::Integer,
    full_slopes_schema::AbstractString,
    selected_slopes_schema::AbstractString,
    lenslet_order_schema::AbstractString,
)
    lenslets = Int(n_lenslets)
    selected_count = Int(selected_lenslet_count)
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann slope-selection n_lenslets must be positive",
    ))
    full_count = Base.checked_mul(lenslets, lenslets)
    0 < selected_count <= full_count || throw(AlgorithmGraphError(
        "Shack–Hartmann selected_lenslet_count must lie in 1:$full_count",
    ))
    schemas = (
        full_slopes_schema,
        selected_slopes_schema,
        lenslet_order_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "Shack–Hartmann slope-selection schemas must not be empty",
    ))
    return ShackHartmannSlopeSelectionNodeConfig(
        lenslets,
        selected_count,
        String(full_slopes_schema),
        String(selected_slopes_schema),
        String(lenslet_order_schema),
    )
end

"""
    shack_hartmann_slope_selection_node(name; ...)

Declare one complete-frame slope selector. The node converts the full AOS
`[axis 1; axis 2]` slope vector into pair-interleaved `[axis 1, axis 2]`
pairs using an explicit startup lenslet order. Instrument ROI-file parsing
remains an application concern.
"""
function shack_hartmann_slope_selection_node(
    name::Symbol;
    n_lenslets::Integer,
    selected_lenslet_count::Integer,
    full_slopes_schema::AbstractString,
    selected_slopes_schema::AbstractString,
    lenslet_order_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ShackHartmannSlopeSelectionNodeConfig(
        n_lenslets,
        selected_lenslet_count,
        full_slopes_schema,
        selected_slopes_schema,
        lenslet_order_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannSlopeSelectionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannSlopeSelectionNode{T}},
    config::ShackHartmannSlopeSelectionNodeConfig,
) where {T}
    full_count = config.n_lenslets * config.n_lenslets
    selected_count = config.selected_lenslet_count
    return (
        graph_port_contract(
            :full_slopes,
            :input,
            :data,
            T,
            (2 * full_count,),
            config.full_slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :selected_slopes,
            :output,
            :data,
            T,
            (2 * selected_count,),
            config.selected_slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :lenslet_order,
            :input,
            :parameter,
            UInt32,
            (selected_count,),
            config.lenslet_order_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannSlopeSelectionOwner{P,I,O}
    plan::P
    input::I
    output::O
end

function prepare_graph_node(
    ::Type{ShackHartmannSlopeSelectionNode{T}},
    config::ShackHartmannSlopeSelectionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:full_slopes,)},
    outputs::NamedTuple{(:selected_slopes,)},
    parameters::NamedTuple{(:lenslet_order,)},
    target,
) where {T}
    plan = try
        ShackHartmannSlopeSelectionPlan(
            config.n_lenslets,
            parameters.lenslet_order,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    selected_lenslet_count(plan) == config.selected_lenslet_count || throw(
        AlgorithmGraphError(
            "Shack–Hartmann lenslet order does not match selected_lenslet_count",
        ),
    )
    return _ShackHartmannSlopeSelectionOwner(
        plan,
        inputs.full_slopes,
        outputs.selected_slopes,
    )
end

@inline function step_graph_node!(owner::_ShackHartmannSlopeSelectionOwner)
    select_shack_hartmann_slopes!(
        owner.output,
        owner.plan,
        owner.input,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_ShackHartmannSlopeSelectionOwner)
    fill!(owner.output, zero(eltype(owner.output)))
    return nothing
end

struct CCDDetectorAcquisitionNode{T<:AbstractFloat} end

"""Construction values for one complete-frame single-read CCD acquisition."""
struct CCDDetectorAcquisitionNodeConfig{T<:AbstractFloat,N<:NoiseModel}
    rows::Int
    columns::Int
    binning::Int
    pixel_scale_arcsec::T
    wavelength_m::T
    exposure_duration_s::T
    quantum_efficiency::T
    noise::N
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

@inline function _ccd_detector_noise(
    ::Type{T},
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
) where {T<:AbstractFloat}
    sigma = T(readout_noise_e)
    isfinite(sigma) && sigma >= zero(T) || throw(AlgorithmGraphError(
        "CCD readout_noise_e must be finite and nonnegative",
    ))
    if readout_noise
        sigma > zero(T) || throw(AlgorithmGraphError(
            "CCD readout_noise_e must be positive when readout noise is enabled",
        ))
        return photon_noise ? NoisePhotonReadout{T}(sigma) :
            NoiseReadout{T}(sigma)
    end
    iszero(sigma) || throw(AlgorithmGraphError(
        "CCD readout_noise_e must be zero when readout noise is disabled",
    ))
    return photon_noise ? NoisePhoton() : NoiseNone()
end

function _ccd_detector_acquisition_config(
    ::Type{T};
    rows::Integer,
    columns::Integer,
    binning::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
) where {T<:AbstractFloat}
    n_rows = Int(rows)
    n_columns = Int(columns)
    bin = Int(binning)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "CCD photon-rate dimensions must be positive",
    ))
    bin > 0 || throw(AlgorithmGraphError(
        "CCD binning must be positive",
    ))
    n_rows % bin == 0 && n_columns % bin == 0 || throw(
        AlgorithmGraphError(
            "CCD photon-rate dimensions must be divisible by binning",
        ))

    pixel_scale = T(pixel_scale_arcsec)
    wavelength = T(wavelength_m)
    exposure = T(exposure_duration_s)
    qe = T(quantum_efficiency)
    isfinite(pixel_scale) && pixel_scale > zero(T) || throw(
        AlgorithmGraphError(
            "CCD pixel_scale_arcsec must be finite and positive",
        ))
    isfinite(wavelength) && wavelength > zero(T) || throw(
        AlgorithmGraphError(
            "CCD wavelength_m must be finite and positive",
        ))
    isfinite(exposure) && exposure > zero(T) || throw(
        AlgorithmGraphError(
            "CCD exposure_duration_s must be finite and positive",
        ))
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(
        AlgorithmGraphError(
            "CCD quantum_efficiency must lie in [0, 1]",
        ))
    rng_seed >= 0 || throw(AlgorithmGraphError(
        "CCD rng_seed must be nonnegative",
    ))
    rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "CCD rng_seed exceeds UInt64",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "CCD photon_rate_schema must not be empty",
    ))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "CCD frame_schema must not be empty",
    ))

    noise = _ccd_detector_noise(
        T,
        photon_noise,
        readout_noise,
        readout_noise_e,
    )
    return CCDDetectorAcquisitionNodeConfig(
        n_rows,
        n_columns,
        bin,
        pixel_scale,
        wavelength,
        exposure,
        qe,
        noise,
        UInt64(rng_seed),
        String(photon_rate_schema),
        String(frame_schema),
    )
end

"""
    ccd_detector_acquisition_node(name; ...)

Declare one complete-frame, single-read CCD acquisition. The node consumes a
cell-integrated detector-plane photon-rate mosaic and writes a detector frame.
Its RNG is explicit persistent state and is restored by `reset_graph!`.
Partial exposure timing, rolling shutters, frame transfer, and readout
readiness remain Plant operations.
"""
function ccd_detector_acquisition_node(
    name::Symbol;
    rows::Integer,
    columns::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
    binning::Integer=1,
    photon_noise::Bool=true,
    readout_noise::Bool=false,
    readout_noise_e::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _ccd_detector_acquisition_config(
        T;
        rows,
        columns,
        binning,
        pixel_scale_arcsec,
        wavelength_m,
        exposure_duration_s,
        quantum_efficiency,
        photon_noise,
        readout_noise,
        readout_noise_e,
        rng_seed,
        photon_rate_schema,
        frame_schema,
    )
    return algorithm_node(
        name,
        CCDDetectorAcquisitionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{CCDDetectorAcquisitionNode{T}},
    config::CCDDetectorAcquisitionNodeConfig{T},
) where {T}
    return (
        graph_port_contract(
            :photon_rate,
            :input,
            :data,
            T,
            (config.rows, config.columns),
            config.photon_rate_schema,
            :column_major,
        ),
        graph_port_contract(
            :frame,
            :output,
            :data,
            T,
            (div(config.rows, config.binning),
                div(config.columns, config.binning)),
            config.frame_schema,
            :column_major,
        ),
    )
end

struct _CCDDetectorAcquisitionOwner{P,R,O,A}
    prepared::P
    rng::R
    output::O
    latent_buffer::A
    rng_seed::UInt64
end

@inline _arcsec_to_rad(value::T) where {T<:AbstractFloat} =
    value * T(pi / (180 * 3600))

function prepare_graph_node(
    ::Type{CCDDetectorAcquisitionNode{T}},
    config::CCDDetectorAcquisitionNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:photon_rate,)},
    outputs::NamedTuple{(:frame,)},
    ::NamedTuple{()},
    target,
) where {T}
    sampling = _arcsec_to_rad(config.pixel_scale_arcsec)
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        inputs.photon_rate;
        coordinate_domain=AngularCoordinates(),
        sampling=(sampling, sampling),
        spectral=MonochromaticChannel(config.wavelength_m),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    detector = Detector(
        exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency,
        psf_sampling=1,
        binning=config.binning,
        noise=config.noise,
        sensor=CCDSensor(T=T),
        response_model=NullFrameResponse(),
        T=T,
        backend=compute_device_backend(target),
    )
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(
        AlgorithmGraphError(
            "prepared CCD output shape does not match its graph frame port",
        ))
    rng = runtime_rng(config.rng_seed)
    state = detector_acquisition_state(prepared)
    return _CCDDetectorAcquisitionOwner(
        prepared,
        rng,
        outputs.frame,
        state.latent_buffer,
        config.rng_seed,
    )
end

@inline function step_graph_node!(owner::_CCDDetectorAcquisitionOwner)
    frame = capture!(owner.prepared, owner.rng)
    copyto!(owner.output, frame)
    return nothing
end

function reset_graph_node!(owner::_CCDDetectorAcquisitionOwner)
    detector = detector_acquisition_detector(owner.prepared)
    reset_integration!(detector)
    fill!(owner.latent_buffer, zero(eltype(owner.latent_buffer)))
    fill!(output_frame(detector), zero(eltype(output_frame(detector))))
    fill!(owner.output, zero(eltype(owner.output)))
    seed!(owner.rng, owner.rng_seed)
    return nothing
end

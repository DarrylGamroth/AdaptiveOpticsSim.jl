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

module HILReferenceSystems

using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends: HostComputeDevice
using AdaptiveOpticsSim.Optics
import AdaptiveOpticsSim.Optics: pupil_support
using AdaptiveOpticsSim.Optics: TelescopeDefinition, prepare_telescope

export actuator_coordinates
export actuator_count
export graph_path
export influence_width
export prepare_atmospheric_hil_reference_system
export prepare_hil_reference_system
export prepare_hil_reference_science_diagnostics
export update_hil_reference_science_diagnostics!
export open_loop_psf
export closed_loop_psf
export open_loop_on_axis_strehl
export closed_loop_on_axis_strehl
export shack_hartmann_lenslet_order
export shack_hartmann_reference_signal
export shack_hartmann_valid_subapertures
export supported_wavefront_sensors

const _GRAPH_DIRECTORY = normpath(joinpath(@__DIR__, "..", "graphs"))
const _SUPPORTED_WAVEFRONT_SENSORS = (:shack_hartmann, :pyramid)
const _ACTUATOR_AXIS_COUNT = 5
const _ACTUATOR_COUNT = _ACTUATOR_AXIS_COUNT^2
const _ACTUATOR_PITCH = 0.4
const _ADJACENT_ACTUATOR_COUPLING = 0.3
const _SHACK_HARTMANN_LENSLET_COUNT = 8
const _REFERENCE_RESOLUTION = 64
const _REFERENCE_TELESCOPE_DIAMETER_M = 1.0
const _REFERENCE_SCIENCE_WAVELENGTH_M = 750.0e-9
const _REFERENCE_ATMOSPHERE_OPD_SCHEMA =
    "org.adaptiveopticssim.hil-reference.uncompensated-opd-m.f32/1"

mutable struct PreparedHILReferenceScienceDiagnostics{P,I,A,T}
    open_loop_pupil::P
    closed_loop_pupil::P
    open_loop_imaging::I
    closed_loop_imaging::I
    open_loop_psf::A
    closed_loop_psf::A
    open_loop_coherent_real::A
    open_loop_coherent_imag::A
    closed_loop_coherent_real::A
    closed_loop_coherent_imag::A
    diffraction_limited_peak::T
    diffraction_limited_field_power::T
    wavelength_m::T
    open_loop_on_axis_strehl::T
    closed_loop_on_axis_strehl::T
end

@inline open_loop_psf(diagnostics::PreparedHILReferenceScienceDiagnostics) =
    diagnostics.open_loop_psf
@inline closed_loop_psf(diagnostics::PreparedHILReferenceScienceDiagnostics) =
    diagnostics.closed_loop_psf
@inline function open_loop_on_axis_strehl(
    diagnostics::PreparedHILReferenceScienceDiagnostics,
)
    return diagnostics.open_loop_on_axis_strehl
end
@inline function closed_loop_on_axis_strehl(
    diagnostics::PreparedHILReferenceScienceDiagnostics,
)
    return diagnostics.closed_loop_on_axis_strehl
end
@inline pupil_support(diagnostics::PreparedHILReferenceScienceDiagnostics) =
    pupil_support(diagnostics.open_loop_pupil)

function _prepare_hil_reference_telescope()
    definition = TelescopeDefinition(
        resolution=_REFERENCE_RESOLUTION,
        diameter=_REFERENCE_TELESCOPE_DIAMETER_M,
        central_obstruction=0.0,
        pupil_reflectivity=1.0,
        revision=1,
        T=Float32,
    )
    return prepare_telescope(definition, HostComputeDevice())
end

"""
    prepare_hil_reference_science_diagnostics()

Prepare host-resident, diffraction-limited-normalized, 750 nm science imaging
for the maintained HIL reference aperture. The two paths accept an
uncompensated OPD and its corrected residual independently. Each published PSF
is normalized by the exact, unsampled diffraction-limited peak of the matching
aperture. The on-axis Strehl ratio is evaluated independently from the
coherent pupil sum, so it does not depend on whether the focal-plane sampling
places the optical axis on a pixel. It includes image motion rather than
recentering the PSF.
"""
function prepare_hil_reference_science_diagnostics()
    telescope = _prepare_hil_reference_telescope()
    source = Source(
        band=:custom,
        wavelength=_REFERENCE_SCIENCE_WAVELENGTH_M,
        photon_irradiance=1.0,
        T=Float32,
    )
    open_loop_pupil = PupilFunction(telescope)
    closed_loop_pupil = PupilFunction(telescope)
    open_loop_imaging = prepare_direct_imaging(
        open_loop_pupil,
        source;
        zero_padding=2,
    )
    closed_loop_imaging = prepare_direct_imaging(
        closed_loop_pupil,
        source;
        zero_padding=2,
    )
    open_loop_product = form_direct_image!(open_loop_imaging)
    closed_loop_product = form_direct_image!(closed_loop_imaging)
    open_loop_values = intensity_values(open_loop_product)
    amplitude = pupil_amplitude(open_loop_pupil)
    amplitude_power_values = similar(amplitude)
    @. amplitude_power_values = abs2(amplitude)
    amplitude_power = sum(amplitude_power_values)
    diffraction_limited_field_power = abs2(sum(amplitude))
    diffraction_limited_peak = sum(open_loop_values) *
        diffraction_limited_field_power /
        (length(open_loop_values) * amplitude_power)
    isfinite(diffraction_limited_peak) &&
        diffraction_limited_peak > zero(diffraction_limited_peak) || error(
        "the HIL reference diffraction-limited peak must be finite and positive",
    )
    open_loop_psf_values = similar(open_loop_values)
    closed_loop_psf_values = similar(intensity_values(closed_loop_product))
    open_loop_coherent_real = similar(amplitude)
    open_loop_coherent_imag = similar(amplitude)
    closed_loop_coherent_real = similar(amplitude)
    closed_loop_coherent_imag = similar(amplitude)
    copyto!(open_loop_psf_values, open_loop_values)
    copyto!(closed_loop_psf_values, intensity_values(closed_loop_product))
    open_loop_psf_values ./= diffraction_limited_peak
    closed_loop_psf_values ./= diffraction_limited_peak
    one_strehl = one(diffraction_limited_peak)
    return PreparedHILReferenceScienceDiagnostics(
        open_loop_pupil,
        closed_loop_pupil,
        open_loop_imaging,
        closed_loop_imaging,
        open_loop_psf_values,
        closed_loop_psf_values,
        open_loop_coherent_real,
        open_loop_coherent_imag,
        closed_loop_coherent_real,
        closed_loop_coherent_imag,
        diffraction_limited_peak,
        diffraction_limited_field_power,
        Float32(_REFERENCE_SCIENCE_WAVELENGTH_M),
        one_strehl,
        one_strehl,
    )
end

"""
    update_hil_reference_science_diagnostics!(diagnostics,
        uncompensated_opd, residual_opd)

Form the open-loop and corrected science PSFs from complete pupil OPD maps in
metres. This mutates only the prepared diagnostic owner and returns it.
"""
function update_hil_reference_science_diagnostics!(
    diagnostics::PreparedHILReferenceScienceDiagnostics,
    uncompensated_opd::AbstractMatrix,
    residual_opd::AbstractMatrix,
)
    apply_opd!(diagnostics.open_loop_pupil, uncompensated_opd)
    apply_opd!(diagnostics.closed_loop_pupil, residual_opd)
    open_loop_product = form_direct_image!(diagnostics.open_loop_imaging)
    closed_loop_product = form_direct_image!(diagnostics.closed_loop_imaging)
    open_loop_values = intensity_values(open_loop_product)
    closed_loop_values = intensity_values(closed_loop_product)
    peak = diagnostics.diffraction_limited_peak
    @. diagnostics.open_loop_psf = open_loop_values / peak
    @. diagnostics.closed_loop_psf = closed_loop_values / peak
    amplitude = pupil_amplitude(diagnostics.open_loop_pupil)
    phase_per_opd = Float32(2 * pi) / diagnostics.wavelength_m
    @. diagnostics.open_loop_coherent_real =
        amplitude * cos(phase_per_opd * uncompensated_opd)
    @. diagnostics.open_loop_coherent_imag =
        amplitude * sin(phase_per_opd * uncompensated_opd)
    @. diagnostics.closed_loop_coherent_real =
        amplitude * cos(phase_per_opd * residual_opd)
    @. diagnostics.closed_loop_coherent_imag =
        amplitude * sin(phase_per_opd * residual_opd)
    open_loop_real = sum(diagnostics.open_loop_coherent_real)
    open_loop_imag = sum(diagnostics.open_loop_coherent_imag)
    closed_loop_real = sum(diagnostics.closed_loop_coherent_real)
    closed_loop_imag = sum(diagnostics.closed_loop_coherent_imag)
    reference_power = diagnostics.diffraction_limited_field_power
    diagnostics.open_loop_on_axis_strehl =
        (open_loop_real^2 + open_loop_imag^2) / reference_power
    diagnostics.closed_loop_on_axis_strehl =
        (closed_loop_real^2 + closed_loop_imag^2) / reference_power
    return diagnostics
end

"""Return the WFS families with maintained, fully declared HIL reference systems."""
supported_wavefront_sensors() = _SUPPORTED_WAVEFRONT_SENSORS

"""Return the number of commands in the analytic reference deformable mirror."""
actuator_count() = _ACTUATOR_COUNT

"""
    influence_width([T=Float32])

Return the Gaussian influence width in normalized pupil coordinates. The
reference deformable mirror declares an adjacent-actuator pitch of `0.4` and a
unit-command mechanical coupling of `0.3`, so the width is fixed by
`coupling = exp(-pitch^2 / (2 * width^2))`.
"""
function influence_width(::Type{T}=Float32) where {T<:AbstractFloat}
    pitch = T(_ACTUATOR_PITCH)
    coupling = T(_ADJACENT_ACTUATOR_COUPLING)
    return pitch / sqrt(-T(2) * log(coupling))
end

"""
    actuator_coordinates([T=Float32])

Return the fully declared 2-by-25 normalized-pupil coordinates of the analytic
reference deformable mirror. Command order is column-major on the regular
5-by-5 grid: normalized-pupil x is the fast axis and y is the slow axis.
"""
function actuator_coordinates(
    ::Type{T}=Float32,
) where {T<:AbstractFloat}
    coordinates = Matrix{T}(undef, 2, _ACTUATOR_COUNT)
    centre = (_ACTUATOR_AXIS_COUNT + 1) ÷ 2
    pitch = T(_ACTUATOR_PITCH)
    command_index = 1
    for y_index in 1:_ACTUATOR_AXIS_COUNT
        for x_index in 1:_ACTUATOR_AXIS_COUNT
            coordinates[1, command_index] = T(x_index - centre) * pitch
            coordinates[2, command_index] = T(y_index - centre) * pitch
            command_index += 1
        end
    end
    return coordinates
end

"""
    shack_hartmann_valid_subapertures()

Return the exact 8-by-8 valid-lenslet mask for the SHWFS HIL reference system.
A lenslet is valid when its normalized-pupil centre lies on or inside the unit
pupil. The mask is a declared reference-system rule, not a measured ROI map.
"""
function shack_hartmann_valid_subapertures()
    count = _SHACK_HARTMANN_LENSLET_COUNT
    radius = count / 2
    centre = (count + 1) / 2
    mask = Matrix{Bool}(undef, count, count)
    for lenslet_y in axes(mask, 2), lenslet_x in axes(mask, 1)
        x = (lenslet_x - centre) / radius
        y = (lenslet_y - centre) / radius
        mask[lenslet_x, lenslet_y] = x * x + y * y <= 1
    end
    return mask
end

"""Return the zero centroid reference for all 64 declared SHWFS lenslets."""
shack_hartmann_reference_signal(::Type{T}=Float32) where {T<:AbstractFloat} =
    zeros(T, _SHACK_HARTMANN_LENSLET_COUNT^2, 2)

"""
    shack_hartmann_lenslet_order([T=UInt32])

Return valid lenslets in canonical Julia column-major lenslet order. The
SHWFS slope-selection node publishes interleaved axis-1/axis-2 pairs in this
order.
"""
function shack_hartmann_lenslet_order(
    ::Type{T}=UInt32,
) where {T<:Integer}
    mask = shack_hartmann_valid_subapertures()
    order = Vector{T}(undef, count(mask))
    output_index = 1
    for lenslet_index in eachindex(mask)
        mask[lenslet_index] || continue
        order[output_index] = T(lenslet_index)
        output_index += 1
    end
    return order
end

@inline _graph_filename(::Val{:shack_hartmann}) =
    "shack_hartmann_hil_reference.toml"
@inline _graph_filename(::Val{:pyramid}) = "pyramid_hil_reference.toml"

"""Return the maintained graph file for one fully declared HIL reference system."""
function graph_path(wavefront_sensor::Symbol)
    wavefront_sensor in _SUPPORTED_WAVEFRONT_SENSORS || throw(ArgumentError(
        "unsupported HIL reference WFS '$wavefront_sensor'; expected one of " *
        "$(join(_SUPPORTED_WAVEFRONT_SENSORS, ", "))",
    ))
    return joinpath(_GRAPH_DIRECTORY, _graph_filename(Val(wavefront_sensor)))
end

function _bindings(::Val{:shack_hartmann}, ::Type{T}) where {T<:AbstractFloat}
    return (
        dm_command=zeros(T, _ACTUATOR_COUNT),
        uncompensated_opd=zeros(T, 64, 64),
        dm_actuator_coordinates=actuator_coordinates(T),
        valid_subapertures=shack_hartmann_valid_subapertures(),
        reference_signal=shack_hartmann_reference_signal(T),
        lenslet_order=shack_hartmann_lenslet_order(),
    )
end

function _bindings(::Val{:pyramid}, ::Type{T}) where {T<:AbstractFloat}
    return (
        dm_command=zeros(T, _ACTUATOR_COUNT),
        uncompensated_opd=zeros(T, 64, 64),
        dm_actuator_coordinates=actuator_coordinates(T),
    )
end

"""
    prepare_hil_reference_system(wavefront_sensor; target=HostComputeDevice())

Prepare one fully declared, noiseless HIL reference system and its complete
frame/command lockstep boundary. The returned `uncompensated_opd` is the exact
caller-owned graph input. Mutate it only between completed HIL exchanges.
"""
function prepare_hil_reference_system(
    wavefront_sensor::Symbol;
    target=HostComputeDevice(),
)
    wavefront_sensor in _SUPPORTED_WAVEFRONT_SENSORS ||
        graph_path(wavefront_sensor)
    bindings = _bindings(Val(wavefront_sensor), Float32)
    definition = load_algorithm_graph(
        graph_path(wavefront_sensor);
        bindings,
    )
    graph = prepare_algorithm_graph(definition; target)
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:dm_command,
        frame_output=:wfs_frame,
    )
    return (
        graph,
        boundary,
        uncompensated_opd=graph_input(graph, Val(:uncompensated_opd)),
    )
end

function _reference_atmosphere_node(; atmosphere_step::Real, rng_seed::Integer)
    return multilayer_atmosphere_opd_node(
        :atmosphere;
        resolution=_REFERENCE_RESOLUTION,
        telescope_diameter_m=_REFERENCE_TELESCOPE_DIAMETER_M,
        central_obstruction_ratio=0.0,
        pupil_reflectivity=1.0,
        aperture_revision=1,
        r0=0.30,
        reference_wavelength_m=500e-9,
        L0=30.0,
        fractional_cn2=(0.55, 0.20, 0.15, 0.10),
        wind_speed=(5.0, 7.0, 10.0, 12.0),
        wind_direction_deg=(0.0, 75.0, 170.0, 260.0),
        altitude=(0.0, 2_000.0, 7_000.0, 12_000.0),
        layer_ids=(
            :reference_ground,
            :reference_2km,
            :reference_7km,
            :reference_12km,
        ),
        atmosphere_step,
        rng_seed,
        atmosphere_opd_schema=_REFERENCE_ATMOSPHERE_OPD_SCHEMA,
        T=Float32,
    )
end

"""
    prepare_atmospheric_hil_reference_system(wavefront_sensor;
        atmosphere_step=1e-3, rng_seed=1, target=HostComputeDevice())

Prepare the maintained reference WFS and deformable mirror with a graph-owned,
deterministically evolving four-layer atmosphere. The returned boundary owns
the only command input and publishes complete detector frames in lockstep.

Calibrate an external RTC against [`prepare_hil_reference_system`](@ref), whose
uncompensated OPD input remains flat, before running this atmosphere-backed
graph. Both definitions retain identical DM, WFS, and detector contracts.
"""
function prepare_atmospheric_hil_reference_system(
    wavefront_sensor::Symbol;
    atmosphere_step::Real=1.0e-3,
    rng_seed::Integer=1,
    target=HostComputeDevice(),
)
    wavefront_sensor in _SUPPORTED_WAVEFRONT_SENSORS ||
        graph_path(wavefront_sensor)
    bindings = _bindings(Val(wavefront_sensor), Float32)
    base = load_algorithm_graph(
        graph_path(wavefront_sensor);
        bindings,
    )
    length(base.inputs) == 2 || error(
        "the maintained HIL reference graph must declare command and OPD inputs",
    )
    atmosphere = _reference_atmosphere_node(; atmosphere_step, rng_seed)
    definition = algorithm_graph(
        (atmosphere, base.nodes...);
        name=Symbol(wavefront_sensor, "_atmosphere_hil_reference"),
        inputs=(first(base.inputs),),
        outputs=(
            base.outputs...,
            graph_output(:atmosphere_opd, :atmosphere => :atmosphere_opd),
        ),
        links=(
            link(
                :atmosphere => :atmosphere_opd,
                :pupil_opd_composition => :uncompensated_opd,
            ),
            base.links...,
        ),
        delayed_links=base.delayed_links,
        parameters=base.parameters,
    )
    graph = prepare_algorithm_graph(definition; target)
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:dm_command,
        frame_output=:wfs_frame,
    )
    return (; graph, boundary)
end

end # module HILReferenceSystems

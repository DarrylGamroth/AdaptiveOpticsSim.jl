module REVOLTClassicHIL

using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends: HostComputeDevice
using ..HILReferenceSystems: prepare_hil_science_diagnostics
using ..REVOLTHILGraphs
using ..REVOLTHSDM277

export command_count
export prepare_calibration_system
export prepare_hil_system
export prepare_science_diagnostics
export valid_subapertures

const _COMMAND_COUNT = 277
const _PUPIL_RESOLUTION = 240
const _TELESCOPE_DIAMETER_M = 1.22
const _CENTRAL_OBSTRUCTION_RATIO = 0.20491803278688525
const _SCIENCE_WAVELENGTH_M = 750.0e-9
const _LENSLET_ROW_WIDTHS = (
    4,
    8,
    10,
    12,
    14,
    14,
    16,
    16,
    16,
    16,
    14,
    14,
    12,
    10,
    8,
    4,
)

"""Return the number of physical HSDM277 command elements."""
command_count() = _COMMAND_COUNT

"""
    valid_subapertures()

Return the exact 16-by-16 REVOLT Classic valid-subaperture mask. The 188 valid
lenslets reproduce `validSubapMask.fits` from the maintained REVOLT on-sky
configuration (raw SHA-256
`720bf7aa82d098e1bc5426b3364dd509adc609617b23127351d4d750b0ff596e`).
The source artifact encodes valid and invalid entries as `1` and `-1`; this
function exposes only the corresponding Boolean selection.
"""
function valid_subapertures()
    mask = falses(16, 16)
    for (row, width) in pairs(_LENSLET_ROW_WIDTHS)
        first_column = (size(mask, 2) - width) ÷ 2 + 1
        last_column = first_column + width - 1
        @views fill!(mask[row, first_column:last_column], true)
    end
    return mask
end

function _bindings(profile::Symbol)
    pdm_command = zeros(Float32, _COMMAND_COUNT)
    if profile === :coordinate_gaussian
        return (;
            pdm_command,
            pdm_actuator_coordinates=
                REVOLTHSDM277.actuator_coordinates(Float32),
        )
    elseif profile === :grid_gaussian
        return (;
            pdm_command,
            pdm_actuator_grid_indices=
                REVOLTHSDM277.actuator_grid_indices(Int32),
        )
    end
    throw(ArgumentError(
        "unsupported REVOLT Classic graph profile '$profile'; expected " *
        "one of $(join(REVOLTHILGraphs.supported_profiles(), ", "))",
    ))
end

function _definition(profile::Symbol)
    bindings = _bindings(profile)
    return load_algorithm_graph(
        REVOLTHILGraphs.graph_path(:classic, profile);
        bindings,
    )
end

"""
    prepare_hil_system(; profile=:grid_gaussian, target=HostComputeDevice())

Prepare the maintained, atmosphere-backed REVOLT Classic detector-frame graph
and its 277-command lockstep HIL boundary. `profile=:grid_gaussian` selects the
regular-grid separable evaluation of the provisional Gaussian HSDM277 surface;
it retains the 240-sample pupil, 16-by-16 diffractive Shack–Hartmann optics,
noisy 352-by-352 C-BLUE One IMX425 frame, and exact physical command order.
"""
function prepare_hil_system(;
    profile::Symbol=:grid_gaussian,
    target=HostComputeDevice(),
)
    graph = prepare_algorithm_graph(_definition(profile); target)
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
    )
    return (; graph, boundary)
end

function _calibration_detector(production_detector)
    config = production_detector.config
    return cmos_detector_acquisition_node(
        :detector;
        rows=config.rows,
        columns=config.columns,
        binning=config.binning,
        pixel_scale_arcsec=config.pixel_scale_arcsec,
        wavelength_m=config.wavelength_m,
        exposure_duration_s=config.exposure_duration_s,
        quantum_efficiency=config.quantum_efficiency,
        gain=config.gain,
        dark_current_e_per_pixel_s=0,
        bits=config.bits,
        full_well_e=config.full_well_e,
        photon_noise=false,
        readout_noise=false,
        readout_noise_e=0,
        column_readout_noise_e=0,
        row_readout_noise_e=0,
        rng_seed=config.rng_seed,
        photon_rate_schema=config.photon_rate_schema,
        frame_schema=config.frame_schema,
        T=Float32,
    )
end

"""
    prepare_calibration_system(; profile=:grid_gaussian,
        target=HostComputeDevice())

Prepare a flat, noiseless REVOLT Classic calibration graph. It retains the
selected HSDM277 model and production Shack–Hartmann/CMOS geometry while
replacing the evolving atmosphere with an explicit zero uncompensated pupil
OPD and disabling detector noise. It is intended only for a simulation-local
interaction matrix; it is not an instrument calibration.
"""
function prepare_calibration_system(;
    profile::Symbol=:grid_gaussian,
    target=HostComputeDevice(),
)
    production = _definition(profile)
    length(production.nodes) == 5 || error(
        "the maintained REVOLT Classic HIL graph must contain five nodes",
    )
    pdm = production.nodes[2]
    composition = production.nodes[3]
    shwfs = production.nodes[4]
    detector = _calibration_detector(production.nodes[5])
    uncompensated_opd = zeros(Float32, _PUPIL_RESOLUTION, _PUPIL_RESOLUTION)
    definition = algorithm_graph(
        (pdm, composition, shwfs, detector);
        name=:revolt_classic_hil_calibration,
        inputs=(
            first(production.inputs),
            graph_input(
                :uncompensated_opd,
                :pupil_opd_composition => :uncompensated_opd,
                uncompensated_opd,
            ),
        ),
        outputs=(
            graph_output(:pdm_surface_opd, :pdm => :surface_opd),
            graph_output(:pupil_opd, :pupil_opd_composition => :pupil_opd),
            graph_output(:shwfs_photon_rate, :shwfs => :photon_rate),
            graph_output(:shwfs_frame, :detector => :frame),
        ),
        links=(
            link(:pdm => :surface_opd,
                :pupil_opd_composition => :surface_opd),
            link(:pupil_opd_composition => :pupil_opd, :shwfs => :opd),
            link(:shwfs => :photon_rate, :detector => :photon_rate),
        ),
        parameters=production.parameters,
    )
    graph = prepare_algorithm_graph(definition; target)
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
    )
    return (; graph, boundary, uncompensated_opd)
end

"""Prepare 750 nm science diagnostics for the REVOLT Classic aperture."""
function prepare_science_diagnostics()
    return prepare_hil_science_diagnostics(;
        resolution=_PUPIL_RESOLUTION,
        telescope_diameter_m=_TELESCOPE_DIAMETER_M,
        central_obstruction_ratio=_CENTRAL_OBSTRUCTION_RATIO,
        wavelength_m=_SCIENCE_WAVELENGTH_M,
    )
end

end # module REVOLTClassicHIL

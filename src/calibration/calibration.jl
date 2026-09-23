"""
    Calibration

Plant-side owner of simulated response acquisition, physical calibration
observables, and accepted-product materialization for plant-side runtime use.
Reusable numerical calibration methods and products are owned by
AdaptiveOpticsCalibration.
"""
module Calibration

import AdaptiveOpticsCalibration
import AdaptiveOpticsCalibration.ModalBases: KarhunenLoeveBasis
using KernelAbstractions
using LinearAlgebra
using SparseArrays

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration,
    NumericalConditionError,
    UnsupportedAlgorithm,
    fftfreq!,
    fftshift2d!

import ..Backends:
    AbstractArrayBackend,
    AcceleratorStyle,
    CPUBackend,
    ExecutionStyle,
    ScalarCPUStyle,
    _resolve_array_backend,
    backend_fill,
    backend_sum_value,
    compute_device,
    execute_fft_plan!,
    execution_style,
    gpu_backend_name,
    launch_kernel!,
    plan_fft_backend!,
    plan_ifft_backend!

using ..Optics
import ..Optics:
    AbstractSource,
    actuator_coordinates,
    anamorphosis_angle_deg,
    misregistration_component,
    rotation_deg,
    sampled_influence_matrix,
    supports_dm_misregistration_identification,
    topology,
    topology_axis_count

using ..Atmospheres
import ..Atmospheres: phase_spectrum
using ..Detectors
import ..Detectors:
    detector_noise_symbol,
    detector_readout_sigma

import ..WavefrontSensors
using ..WavefrontSensors
import ..WavefrontSensors: AbstractWFS

include("build_backends.jl")
include("reconstructor_products.jl")
include("modal_basis.jl")
include("modal_opd_expansion.jl")
include("interaction_matrix.jl")
include("control_matrix.jl")
include("fitting_error.jl")
include("gain_sensing_camera.jl")
include("misregistration_identification.jl")
include("ad_sensitivities.jl")
include("api.jl")

end # module Calibration

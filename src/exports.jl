# Curated package API.
#
# Export only names used routinely after `using AdaptiveOpticsSim`. Declare
# stable advanced and developer-facing seams `public` so callers can use
# `AdaptiveOpticsSim.name` or import them explicitly without adding them to the
# ordinary user namespace.

export AdaptiveOpticsSimError, InvalidConfiguration, DimensionMismatchError
export UnsupportedAlgorithm, NumericalConditionError
export Backends, Optics, Detectors, Atmospheres, WavefrontSensors, Calibration, Plant

export FidelityProfile, ScientificProfile, FastProfile, default_fidelity_profile
export runtime_rng, deterministic_reference_rng

export NullReconstructor, ModalReconstructor, FactorizedReconstructor, MappedReconstructor
export ControlledReconstructor
export reconstruct!, reconstruct
export DiscreteIntegratorController, VectorDelayLine, shift_delay!

export AbstractExecutionPolicy
export SequentialExecution, ThreadedExecution, BackendStreamExecution
export DeterministicExecution, AcceleratedKernelsExecution, DaggerExecution
export SimulationEnsemble
export runtime_timing

public controller_output, reset_controller!, supports_controller_reset
public run_ensemble!, ensemble_members, execution_policy
public ensemble_ownership_roots, init_ensemble_scheduler, execute_ensemble!
public init_execution_state

export TomographyAtmosphereParams, LGSAsterismParams, LGSWFSParams
export TomographyParams, TomographyDMParams
export ModelBasedTomography, InteractionMatrixTomography
export SimulationSlopes, InterleavedSlopes, InvertedSlopes
export build_reconstructor, assemble_reconstructor_and_fitting
export zenith_angle_deg, wind_direction_deg
export reconstruct_wavefront_map
export dm_commands

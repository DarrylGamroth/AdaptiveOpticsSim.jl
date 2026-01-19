module AdaptiveOptics

__precompile__(true)

using FFTW
using LinearAlgebra
using Logging
using Random
using Statistics

"""
AdaptiveOptics.jl

Julia adaptive optics simulation toolkit (in development).
"""
const PROJECT_STATUS = :in_development

include("Core/errors.jl")
include("Core/types.jl")
include("Core/utils.jl")
include("Core/workspace.jl")
include("Core/parallel.jl")

include("Optics/telescope.jl")
include("Optics/source.jl")
include("Optics/psf.jl")
include("Optics/zernike.jl")

export PROJECT_STATUS
export AdaptiveOpticsError, InvalidConfiguration, DimensionMismatchError
export Workspace, ensure_psf_buffers!, backend_type
export Telescope, TelescopeParams, TelescopeState, generate_pupil!, reset_opd!, apply_opd!
export Source, SourceParams, wavelength
export compute_psf!, ensure_psf_state!
export ZernikeBasis, compute_zernike!, noll_to_nm
export AbstractOpticalElement, AbstractSource, AbstractAtmosphere, AbstractWFS
export AbstractDetector, AbstractDeformableMirror, SensingMode, Diffractive, Geometric
export ParallelConfig, with_parallel_config

end # module AdaptiveOptics

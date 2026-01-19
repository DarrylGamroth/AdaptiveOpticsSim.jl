abstract type AbstractOpticalElement end
abstract type AbstractSource <: AbstractOpticalElement end
abstract type AbstractAtmosphere <: AbstractOpticalElement end
abstract type AbstractWFS <: AbstractOpticalElement end
abstract type AbstractDetector <: AbstractOpticalElement end
abstract type AbstractDeformableMirror <: AbstractOpticalElement end

abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric <: SensingMode end

sensing_mode(::AbstractWFS) = Diffractive()

abstract type NoiseModel end
struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end
struct NoiseReadout <: NoiseModel end
struct NoisePhotonReadout <: NoiseModel end

abstract type DMApplyMode end
struct DMAdditive <: DMApplyMode end
struct DMReplace <: DMApplyMode end

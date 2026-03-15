abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric <: SensingMode end

abstract type WFSNormalization end
struct MeanValidFluxNormalization <: WFSNormalization end
struct IncidenceFluxNormalization <: WFSNormalization end

sensing_mode(::AbstractWFS) = Diffractive()

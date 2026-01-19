abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric <: SensingMode end

sensing_mode(::AbstractWFS) = Diffractive()

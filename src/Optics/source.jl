const BAND_WAVELENGTHS = Dict(
    :U => 0.365e-6,
    :B => 0.445e-6,
    :V => 0.551e-6,
    :R => 0.658e-6,
    :I => 0.806e-6,
    :J => 1.22e-6,
    :H => 1.63e-6,
    :K => 2.19e-6,
)

struct SourceParams{T<:AbstractFloat}
    band::Symbol
    magnitude::T
    coordinates::NTuple{2,T}
    wavelength::T
end

struct Source{P<:SourceParams} <: AbstractSource
    params::P
end

function Source(; band::Symbol=:I, magnitude::Real=0.0, coordinates=(0.0, 0.0), wavelength=nothing, T::Type{<:AbstractFloat}=Float64)
    if wavelength === nothing
        if haskey(BAND_WAVELENGTHS, band)
            wavelength = BAND_WAVELENGTHS[band]
        else
            throw(InvalidConfiguration("Unknown band $(band); provide wavelength."))
        end
    end
    coords = (T(coordinates[1]), T(coordinates[2]))
    params = SourceParams{T}(band, T(magnitude), coords, T(wavelength))
    return Source(params)
end

wavelength(src::Source) = src.params.wavelength

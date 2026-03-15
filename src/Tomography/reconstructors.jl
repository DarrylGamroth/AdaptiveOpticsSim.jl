using LinearAlgebra

abstract type AbstractTomographyMethod end

struct ModelBasedTomography <: AbstractTomographyMethod end

struct InteractionMatrixTomography <: AbstractTomographyMethod end

struct TomographicReconstructor{
    Method<:AbstractTomographyMethod,
    T<:AbstractFloat,
    R<:AbstractMatrix{T},
    G<:AbstractMatrix{Bool},
    AP<:TomographyAtmosphereParams,
    LP<:LGSAsterismParams,
    WP<:LGSWFSParams,
    TP<:TomographyParams,
    DP<:TomographyDMParams,
    F,
}
    method::Method
    reconstructor::R
    grid_mask::G
    atmosphere::AP
    asterism::LP
    wfs::WP
    tomography::TP
    dm::DP
    fitting::F
end

function build_reconstructor(
    ::InteractionMatrixTomography,
    interaction_matrix::AbstractMatrix{T},
    grid_mask::AbstractMatrix{Bool},
    atmosphere::TomographyAtmosphereParams,
    asterism::LGSAsterismParams,
    wfs::LGSWFSParams,
    tomography::TomographyParams,
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
    rcond::Real=sqrt(eps(T)),
) where {T<:AbstractFloat}
    size(interaction_matrix, 1) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one row"))
    size(interaction_matrix, 2) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one column"))
    recon = pinv(Matrix(interaction_matrix); rtol=T(rcond))
    return TomographicReconstructor(
        InteractionMatrixTomography(),
        recon,
        Matrix(grid_mask),
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm,
        fitting,
    )
end

function build_reconstructor(
    ::InteractionMatrixTomography,
    imat::InteractionMatrix,
    grid_mask::AbstractMatrix{Bool},
    atmosphere::TomographyAtmosphereParams,
    asterism::LGSAsterismParams,
    wfs::LGSWFSParams,
    tomography::TomographyParams,
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
    rcond::Real=sqrt(eps(eltype(imat.matrix))),
)
    return build_reconstructor(
        InteractionMatrixTomography(),
        imat.matrix,
        grid_mask,
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm;
        fitting=fitting,
        rcond=rcond,
    )
end

function build_reconstructor(
    ::ModelBasedTomography,
    atmosphere::TomographyAtmosphereParams,
    asterism::LGSAsterismParams,
    wfs::LGSWFSParams,
    tomography::TomographyParams,
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
)
    throw(UnsupportedAlgorithm(
        "model-based tomography is not implemented yet; port Cxx/Cox/Cnz/RecStatSA next",
    ))
end

function reconstruct_wavefront!(
    out::AbstractVector{T},
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    size(reconstructor.reconstructor, 1) == length(out) ||
        throw(DimensionMismatchError("output length must match reconstructor row count"))
    size(reconstructor.reconstructor, 2) == length(slopes) ||
        throw(DimensionMismatchError("slopes length must match reconstructor column count"))
    mul!(out, reconstructor.reconstructor, slopes)
    return out
end

function reconstruct_wavefront(
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    out = Vector{T}(undef, size(reconstructor.reconstructor, 1))
    return reconstruct_wavefront!(out, reconstructor, slopes)
end

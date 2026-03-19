using Random
using Statistics

abstract type NCPABasis end
struct KLBasis{M<:KLBasisMethod} <: NCPABasis
    method::M
end
struct ZernikeModalBasis <: NCPABasis end
struct M2CBasis <: NCPABasis end

KLBasis() = KLBasis(KLHHtPSD())

@kernel function combine_basis_kernel!(opd, basis, coeffs, pupil, n_modes::Int)
    I = @index(Global, Cartesian)
    i, j = Tuple(I)
    if i <= size(opd, 1) && j <= size(opd, 2)
        acc = zero(eltype(opd))
        @inbounds for k in 1:n_modes
            acc += coeffs[k] * basis[i, j, k]
        end
        @inbounds opd[i, j] = ifelse(pupil[i, j], acc, zero(acc))
    end
end

struct NCPA{T<:AbstractFloat,A<:AbstractMatrix{T},B,V} <: AbstractOpticalElement
    opd::A
    basis::B
    coeffs::V
end

function NCPA(tel::Telescope, dm::DeformableMirror, atm::AbstractAtmosphere;
    basis::NCPABasis=KLBasis(), coefficients=nothing, f2=nothing, seed::Integer=5,
    M2C::Union{Nothing,AbstractMatrix}=nothing)

    T = eltype(tel.state.opd)
    if f2 === nothing
        if coefficients === nothing
            opd = copy(tel.state.pupil)
            return NCPA{T, typeof(opd), Nothing, Nothing}(opd, nothing, nothing)
        end
        coeffs = T.(coefficients)
        basis_grid = ncpa_basis(basis, tel, dm, atm; n_modes=length(coeffs), M2C=M2C)
        opd = combine_basis(basis_grid, coeffs, tel.state.pupil)
        return NCPA{T, typeof(opd), typeof(basis_grid), typeof(coeffs)}(opd, basis_grid, coeffs)
    end

    length(f2) == 4 ||
        throw(InvalidConfiguration("f2 must be (amplitude, start_mode, end_mode, cutoff)"))
    amplitude, start_mode, end_mode, cutoff = f2
    n_modes = Int(end_mode)
    basis_grid = ncpa_basis(basis, tel, dm, atm; n_modes=n_modes, M2C=M2C)
    rng = MersenneTwister(seed)
    coeffs = zeros(T, n_modes)
    for i in Int(start_mode):Int(end_mode)
        coeffs[i] = T(randn(rng)) / sqrt(T(i) + T(cutoff))
    end
    opd = combine_basis(basis_grid, coeffs, tel.state.pupil)
    σ = std(opd[tel.state.pupil])
    if σ > 0
        opd .*= T(amplitude) / σ
    end
    return NCPA{T, typeof(opd), typeof(basis_grid), typeof(coeffs)}(opd, basis_grid, coeffs)
end

function ncpa_basis(basis::KLBasis{<:KLDMModes}, tel::Telescope, dm::DeformableMirror,
    ::AbstractAtmosphere; n_modes::Int=1, M2C::Union{Nothing,AbstractMatrix}=nothing)
    _, basis_grid = kl_modal_basis(KLDMModes(), dm, tel; n_modes=n_modes)
    return basis_grid
end

function ncpa_basis(basis::KLBasis{<:KLDMModes}, tel::Telescope, dm::DeformableMirror;
    n_modes::Int=1, M2C::Union{Nothing,AbstractMatrix}=nothing)
    _, basis_grid = kl_modal_basis(KLDMModes(), dm, tel; n_modes=n_modes)
    return basis_grid
end

function ncpa_basis(basis::KLBasis{<:KLHHtPSD}, tel::Telescope, dm::DeformableMirror,
    atm::AbstractAtmosphere; n_modes::Int=1, M2C::Union{Nothing,AbstractMatrix}=nothing)
    _, basis_grid = kl_modal_basis(basis.method, dm, tel, atm; n_modes=n_modes)
    return basis_grid
end

function ncpa_basis(::KLBasis{<:KLHHtPSD}, tel::Telescope, dm::DeformableMirror;
    n_modes::Int=1, M2C::Union{Nothing,AbstractMatrix}=nothing)
    throw(InvalidConfiguration("KLBasis(KLHHtPSD()) requires an atmosphere"))
end

function ncpa_basis(::ZernikeModalBasis, tel::Telescope, dm::DeformableMirror; n_modes::Int=1,
    M2C::Union{Nothing,AbstractMatrix}=nothing)
    zb = ZernikeBasis(tel, n_modes)
    compute_zernike!(zb, tel)
    return zb.modes
end

function ncpa_basis(::M2CBasis, tel::Telescope, dm::DeformableMirror; n_modes::Int=1,
    M2C::Union{Nothing,AbstractMatrix}=nothing)
    if M2C === nothing
        throw(InvalidConfiguration("M2C matrix must be provided when basis = M2CBasis()"))
    end
    return basis_from_m2c(dm, tel, M2C)
end

function ncpa_basis(basis::NCPABasis, tel::Telescope, dm::DeformableMirror,
    ::AbstractAtmosphere; n_modes::Int=1, M2C::Union{Nothing,AbstractMatrix}=nothing)
    return ncpa_basis(basis, tel, dm; n_modes=n_modes, M2C=M2C)
end

function combine_basis!(opd::AbstractMatrix{T}, basis::AbstractArray{T,3},
    coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return combine_basis!(execution_style(opd), opd, basis, coeffs, pupil)
end

function combine_basis!(::ScalarCPUStyle, opd::AbstractMatrix{T}, basis::AbstractArray{T,3},
    coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    fill!(opd, zero(T))
    @inbounds for k in 1:n_modes
        @views @. opd += coeffs[k] * basis[:, :, k]
    end
    @. opd *= pupil
    return opd
end

function combine_basis!(style::AcceleratorStyle, opd::AbstractMatrix{T}, basis::AbstractArray{T,3},
    coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    if n_modes == 0
        fill!(opd, zero(T))
        return opd
    end
    launch_kernel!(style, combine_basis_kernel!, opd, basis, coeffs, pupil, n_modes; ndrange=size(opd))
    return opd
end

function combine_basis(basis::AbstractArray{T,3}, coeffs::AbstractVector{T},
    pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n, m = size(basis, 1), size(basis, 2)
    opd = similar(basis, T, n, m)
    return combine_basis!(opd, basis, coeffs, pupil)
end

function apply!(ncpa::NCPA, tel::Telescope, ::DMAdditive)
    if size(ncpa.opd) != size(tel.state.opd)
        throw(DimensionMismatchError("NCPA OPD size does not match telescope resolution"))
    end
    tel.state.opd .+= ncpa.opd
    return tel
end

function apply!(ncpa::NCPA, tel::Telescope, ::DMReplace)
    if size(ncpa.opd) != size(tel.state.opd)
        throw(DimensionMismatchError("NCPA OPD size does not match telescope resolution"))
    end
    tel.state.opd .= ncpa.opd
    return tel
end

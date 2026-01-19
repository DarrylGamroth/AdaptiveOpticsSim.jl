using Random
using Statistics

struct NCPA{T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractArray{T,3},V<:AbstractVector{T}} <: AbstractOpticalElement
    opd::A
    basis::Union{Nothing,B}
    coeffs::Union{Nothing,V}
end

function NCPA(tel::Telescope, dm::DeformableMirror, atm::AbstractAtmosphere;
    modal_basis::Symbol=:KL, coefficients=nothing, f2=nothing, seed::Integer=5,
    M2C::Union{Nothing,AbstractMatrix}=nothing)

    T = eltype(tel.state.opd)
    if f2 === nothing
        if coefficients === nothing
            opd = copy(tel.state.pupil)
            return NCPA{T, typeof(opd), Array{T,3}, Vector{T}}(opd, nothing, nothing)
        end
        coeffs = T.(coefficients)
        basis = ncpa_basis(tel, dm; modal_basis=modal_basis, n_modes=length(coeffs), M2C=M2C)
        opd = combine_basis(basis, coeffs, tel.state.pupil)
        return NCPA{T, typeof(opd), typeof(basis), typeof(coeffs)}(opd, basis, coeffs)
    end

    f2_params = collect(f2)
    if length(f2_params) != 4
        throw(InvalidConfiguration("f2 must be (amplitude, start_mode, end_mode, cutoff)"))
    end
    amplitude, start_mode, end_mode, cutoff = f2_params
    n_modes = Int(end_mode)
    basis = ncpa_basis(tel, dm; modal_basis=modal_basis, n_modes=n_modes, M2C=M2C)
    rng = MersenneTwister(seed)
    coeffs = zeros(T, n_modes)
    for i in Int(start_mode):Int(end_mode)
        coeffs[i] = T(randn(rng)) / sqrt(T(i) + T(cutoff))
    end
    opd = combine_basis(basis, coeffs, tel.state.pupil)
    σ = std(opd[tel.state.pupil])
    if σ > 0
        opd .*= T(amplitude) / σ
    end
    return NCPA{T, typeof(opd), typeof(basis), typeof(coeffs)}(opd, basis, coeffs)
end

function ncpa_basis(tel::Telescope, dm::DeformableMirror; modal_basis::Symbol=:KL, n_modes::Int=1,
    M2C::Union{Nothing,AbstractMatrix}=nothing)
    if modal_basis == :KL
        _, basis = kl_modal_basis(dm, tel; n_modes=n_modes)
        return basis
    elseif modal_basis == :Zernike
        zb = ZernikeBasis(tel, n_modes)
        compute_zernike!(zb, tel)
        return zb.modes
    elseif modal_basis == :M2C
        if M2C === nothing
            throw(InvalidConfiguration("M2C matrix must be provided when modal_basis = :M2C"))
        end
        return basis_from_m2c(dm, tel, M2C)
    end
    throw(InvalidConfiguration("Unsupported modal_basis: $(modal_basis)"))
end

function combine_basis(basis::AbstractArray{T,3}, coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    n, m = size(basis, 1), size(basis, 2)
    opd = zeros(T, n, m)
    @inbounds for k in 1:n_modes
        opd .+= coeffs[k] .* basis[:, :, k]
    end
    opd .*= pupil
    return opd
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

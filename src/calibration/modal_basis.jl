using LinearAlgebra
using Statistics

#
# Modal basis construction
#
# This file builds the command-basis operators used for modal calibration and
# control.
#
# Supported Karhunen-Loeve-like constructions:
# - `KLDMModes`: eigendecomposition of the DM-mode covariance `M' * M`
# - `KLHHtPSD`: PSD-weighted covariance in Fourier space, which better reflects
#   atmospheric statistics when an atmosphere model is available
#
# The resulting `M2C` operator maps modal coefficients to actuator commands,
# while `basis` stores the corresponding OPD modes on the pupil grid.
#
abstract type KLBasisMethod end
struct KLDMModes <: KLBasisMethod end
struct KLHHtPSD <: KLBasisMethod end

"""
    ModalBasis

Bundle the command-space and pupil-space representation of a modal basis.

- `M2C` maps modal coefficients to DM commands
- `basis` stores the corresponding basis vectors in flattened pupil form
- `projector` is an optional inverse/projection operator back into modal space
"""
struct ModalBasis{T<:AbstractFloat,
    M<:AbstractMatrix{T},
    B<:AbstractMatrix{T},
    P<:AbstractMatrix{T}}
    M2C::M
    basis::B
    projector::Union{Nothing,P}
end

@inline modal_to_command(basis::ModalBasis) = basis.M2C
@inline sampled_basis(basis::ModalBasis) = basis.basis
@inline modal_projector(basis::ModalBasis) = basis.projector

function dm_basis(dm::DeformableMirror, tel::Telescope)
    n = tel.params.resolution
    n_modes = size(dm.state.modes, 2)
    return reshape(dm.state.modes, n, n, n_modes)
end

function basis_from_m2c(dm::DeformableMirror, tel::Telescope, M2C::AbstractMatrix)
    n = tel.params.resolution
    basis_mat = dm.state.modes * M2C
    return reshape(basis_mat, n, n, size(M2C, 2))
end

"""
    basis_projector(basis; tol=1e-3, policy=...)

Construct a projector from sampled basis vectors back into modal coefficients.

If the basis is close to diagonal in its Gram matrix, this uses the cheaper
diagonal approximation. Otherwise it falls back to the configured inverse
operator.
"""
function basis_projector(basis::AbstractMatrix{T}; tol::Real=1e-3,
    policy::InversePolicy=default_projector_inverse_policy(T)) where {T<:AbstractFloat}
    cross = transpose(basis) * basis
    diag_vals = diag(cross)
    diag_sum = sum(abs, diag_vals)
    if diag_sum == 0
        projector, _ = inverse_operator(basis, policy)
        return projector
    end
    non_diag_sum = sum(abs, cross) - diag_sum
    criteria = abs(diag_sum - non_diag_sum) / diag_sum
    if criteria <= tol && all(!iszero, diag_vals)
        return Matrix(Diagonal(inv.(diag_vals)) * transpose(basis))
    end
    projector, _ = inverse_operator(basis, policy)
    return projector
end

"""
    kl_modal_basis(method, dm, tel; ...)

Build a KL-style modal basis and its modal-to-command matrix.

The chosen method determines how modal covariance is approximated before the
eigendecomposition that orders modes by decreasing expected variance.
"""
function kl_modal_basis(dm::DeformableMirror, tel::Telescope;
    n_modes::Int=size(dm.state.modes, 2), remove_piston::Bool=true)
    return kl_modal_basis(KLDMModes(), dm, tel; n_modes=n_modes, remove_piston=remove_piston)
end

function kl_modal_basis(::KLDMModes, dm::DeformableMirror, tel::Telescope;
    n_modes::Int=size(dm.state.modes, 2), remove_piston::Bool=true)
    M = dm.state.modes
    C = transpose(M) * M
    evals, evecs = eigen(Symmetric(C))
    order = sortperm(evals; rev=true)
    n_keep = min(n_modes, length(order))
    V = evecs[:, order[1:n_keep]]
    basis_mat = M * V
    if remove_piston
        n = tel.params.resolution
        for i in 1:n_keep
            mode = reshape(view(basis_mat, :, i), n, n)
            mean_mode = mean(mode[tel.state.pupil])
            @views basis_mat[:, i] .-= mean_mode
        end
    end
    basis = reshape(basis_mat, tel.params.resolution, tel.params.resolution, n_keep)
    return V, basis
end

function kl_modal_basis(::KLHHtPSD, dm::DeformableMirror, tel::Telescope, atm::AbstractAtmosphere;
    n_modes::Int=size(dm.state.modes, 2), remove_piston::Bool=true, delta::Union{Nothing,Real}=nothing)
    n = tel.params.resolution
    mode_count = size(dm.state.modes, 2)
    n_keep = min(n_modes, mode_count)
    T = eltype(dm.state.modes)

    delta_val = delta === nothing ? tel.params.diameter / n : delta
    freqs = Vector{T}(undef, n)
    fftfreq!(freqs, n; d=delta_val)
    psd = Matrix{T}(undef, n, n)
    r0, L0 = turbulence_params(atm)
    cst = T(0.023) * T(r0)^(-T(5) / T(3))
    inv_L0 = T(1) / T(L0)
    @inbounds for i in 1:n, j in 1:n
        f = sqrt(freqs[i]^2 + freqs[j]^2)
        k = T(2 * pi) * f
        psd[i, j] = cst * (k^2 + inv_L0^2)^(-T(11) / T(6))
    end
    psd[div(n, 2) + 1, div(n, 2) + 1] = 0

    F = Matrix{Complex{T}}(undef, n * n, mode_count)
    buffer = Matrix{Complex{T}}(undef, n, n)
    fft_plan = plan_fft_backend!(buffer)
    @inbounds for k in 1:mode_count
        mode = reshape(view(dm.state.modes, :, k), n, n)
        @inbounds for i in 1:n, j in 1:n
            buffer[i, j] = complex(mode[i, j], zero(T))
        end
        execute_fft_plan!(buffer, fft_plan)
        @views F[:, k] .= reshape(buffer, :)
    end

    psd_vec = reshape(psd, :)
    weighted = similar(F)
    @inbounds for k in 1:mode_count
        col = view(weighted, :, k)
        fcol = view(F, :, k)
        for idx in 1:length(psd_vec)
            col[idx] = fcol[idx] * psd_vec[idx]
        end
    end

    C = real(transpose(conj(F)) * weighted)
    evals, evecs = eigen(Symmetric(C))
    order = sortperm(evals; rev=true)
    V = evecs[:, order[1:n_keep]]
    basis_mat = dm.state.modes * V
    if remove_piston
        for i in 1:n_keep
            mode = reshape(view(basis_mat, :, i), n, n)
            mean_mode = mean(mode[tel.state.pupil])
            @views basis_mat[:, i] .-= mean_mode
        end
    end
    basis = reshape(basis_mat, n, n, n_keep)
    return V, basis
end

function turbulence_params(atm::KolmogorovAtmosphere)
    return atm.params.r0, atm.params.L0
end

function turbulence_params(atm::MultiLayerAtmosphere)
    return atm.params.r0, atm.params.L0
end

function turbulence_params(atm::AbstractAtmosphere)
    throw(InvalidConfiguration("turbulence_params not defined for $(typeof(atm))"))
end

"""
    modal_basis(dm, tel; ...)

Build the modal command basis used by AO calibration and control.

This returns both the modal-to-command operator and the sampled pupil-space
basis, with an optional projector back into modal coordinates.
"""
function modal_basis(dm::DeformableMirror, tel::Telescope; n_modes::Int=size(dm.state.modes, 2),
    remove_piston::Bool=true, projector::Bool=true, method::KLBasisMethod=KLDMModes(),
    atm::Union{Nothing,AbstractAtmosphere}=nothing)
    M2C, basis = modal_basis_components(method, dm, tel, atm; n_modes=n_modes, remove_piston=remove_piston)
    basis_mat = reshape(basis, :, size(basis, 3))
    proj = projector ? basis_projector(basis_mat) : nothing
    return ModalBasis(M2C, basis_mat, proj)
end

function modal_basis_components(method::KLBasisMethod, dm::DeformableMirror, tel::Telescope, atm;
    n_modes::Int, remove_piston::Bool)
    return kl_modal_basis(method, dm, tel; n_modes=n_modes, remove_piston=remove_piston)
end

function modal_basis_components(method::KLHHtPSD, dm::DeformableMirror, tel::Telescope, ::Nothing;
    n_modes::Int, remove_piston::Bool)
    throw(InvalidConfiguration("KLHHtPSD modal basis requires an atmosphere"))
end

function modal_basis_components(method::KLHHtPSD, dm::DeformableMirror, tel::Telescope, atm::AbstractAtmosphere;
    n_modes::Int, remove_piston::Bool)
    return kl_modal_basis(method, dm, tel, atm; n_modes=n_modes, remove_piston=remove_piston)
end

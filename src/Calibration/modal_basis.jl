using LinearAlgebra
using Statistics

abstract type KLBasisMethod end
struct KLDMModes <: KLBasisMethod end
struct KLHHtPSD <: KLBasisMethod end

struct ModalBasis{T<:AbstractFloat,
    M<:AbstractMatrix{T},
    B<:AbstractMatrix{T},
    P<:AbstractMatrix{T}}
    M2C::M
    basis::B
    projector::Union{Nothing,P}
end

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
        mul!(buffer, fft_plan, buffer)
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

function modal_basis(dm::DeformableMirror, tel::Telescope; n_modes::Int=size(dm.state.modes, 2),
    remove_piston::Bool=true, projector::Bool=true, method::KLBasisMethod=KLDMModes(),
    atm::Union{Nothing,AbstractAtmosphere}=nothing)
    if method isa KLHHtPSD
        if atm === nothing
            throw(InvalidConfiguration("KLHHtPSD modal basis requires an atmosphere"))
        end
        M2C, basis = kl_modal_basis(method, dm, tel, atm; n_modes=n_modes, remove_piston=remove_piston)
    else
        M2C, basis = kl_modal_basis(method, dm, tel; n_modes=n_modes, remove_piston=remove_piston)
    end
    basis_mat = reshape(basis, :, size(basis, 3))
    proj = projector ? basis_projector(basis_mat) : nothing
    return ModalBasis(M2C, basis_mat, proj)
end

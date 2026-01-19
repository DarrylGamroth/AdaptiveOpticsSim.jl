using LinearAlgebra
using Statistics

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

function basis_projector(basis::AbstractMatrix{T}; tol::Real=1e-3) where {T<:AbstractFloat}
    cross = transpose(basis) * basis
    diag_vals = diag(cross)
    diag_sum = sum(abs, diag_vals)
    if diag_sum == 0
        return Matrix(pinv(basis))
    end
    non_diag_sum = sum(abs, cross) - diag_sum
    criteria = abs(diag_sum - non_diag_sum) / diag_sum
    if criteria <= tol && all(!iszero, diag_vals)
        return Matrix(Diagonal(inv.(diag_vals)) * transpose(basis))
    end
    return Matrix(pinv(basis))
end

function kl_modal_basis(dm::DeformableMirror, tel::Telescope; n_modes::Int=size(dm.state.modes, 2), remove_piston::Bool=true)
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

function modal_basis(dm::DeformableMirror, tel::Telescope; n_modes::Int=size(dm.state.modes, 2),
    remove_piston::Bool=true, projector::Bool=true)
    M2C, basis = kl_modal_basis(dm, tel; n_modes=n_modes, remove_piston=remove_piston)
    basis_mat = reshape(basis, :, size(basis, 3))
    proj = projector ? basis_projector(basis_mat) : nothing
    return ModalBasis(M2C, basis_mat, proj)
end

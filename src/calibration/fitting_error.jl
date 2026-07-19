function fitting_error(opd::AbstractMatrix{T}, projector::AbstractMatrix{T}, basis::AbstractMatrix{T}) where {T<:AbstractFloat}
    phi = reshape(opd, :)
    coeffs = projector * phi
    phi_corr = basis * coeffs
    phi_fit = phi .- phi_corr
    n, m = size(opd)
    return reshape(phi_fit, n, m), reshape(phi_corr, n, m), copy(opd)
end

function fitting_error_dm(opd::AbstractMatrix{T}, projector::AbstractMatrix{T},
    pupil::PupilFunction, dm::DeformableMirror,
    M2C::AbstractMatrix{T}) where {T<:AbstractFloat}
    phi = reshape(opd, :)
    coeffs = M2C * (projector * phi)
    dm.state.coefs .= -coeffs
    pupil.opd .= opd
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    opd_fit = copy(pupil.opd)
    opd_corr = dm.state.opd .* pupil.support
    return opd_fit, opd_corr, copy(opd), coeffs
end

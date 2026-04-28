function fitting_error(opd::AbstractMatrix{T}, projector::AbstractMatrix{T}, basis::AbstractMatrix{T}) where {T<:AbstractFloat}
    phi = reshape(opd, :)
    coeffs = projector * phi
    phi_corr = basis * coeffs
    phi_fit = phi .- phi_corr
    n, m = size(opd)
    return reshape(phi_fit, n, m), reshape(phi_corr, n, m), copy(opd)
end

function fitting_error_dm(opd::AbstractMatrix{T}, projector::AbstractMatrix{T},
    tel::Telescope, dm::DeformableMirror, M2C::AbstractMatrix{T}) where {T<:AbstractFloat}
    phi = reshape(opd, :)
    coeffs = M2C * (projector * phi)
    dm.state.coefs .= -coeffs
    tel.state.opd .= opd
    apply!(dm, tel, DMAdditive())
    opd_fit = copy(tel.state.opd)
    opd_corr = dm.state.opd .* tel.state.pupil
    return opd_fit, opd_corr, copy(opd), coeffs
end

# Plant-specific fitting evaluation. Pure explicit-array modal fitting is owned
# by AdaptiveOpticsCalibration.ModalBases.
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

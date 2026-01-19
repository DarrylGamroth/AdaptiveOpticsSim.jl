struct AOCalibration{T<:AbstractFloat}
    basis::Matrix{T}
    projector::Union{Nothing,Matrix{T}}
    M2C::Matrix{T}
    calibration::CalibrationVault{T,Matrix{T}}
end

function ao_calibration(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS;
    n_modes::Int=size(dm.state.modes, 2), amplitude::Real=1e-9,
    projector::Bool=true, basis::Union{Nothing,ModalBasis}=nothing)

    if basis === nothing
        basis = modal_basis(dm, tel; n_modes=n_modes, projector=projector)
    end

    imat = interaction_matrix(dm, wfs, tel, basis.M2C; amplitude=amplitude)
    calib = CalibrationVault(Matrix(imat.matrix))
    return AOCalibration(basis.basis, basis.projector, Matrix(basis.M2C), calib)
end

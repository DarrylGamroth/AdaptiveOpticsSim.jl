struct AOCalibration{T<:AbstractFloat,
    B<:AbstractMatrix{T},
    P<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    C<:CalibrationVault{T}}
    basis::B
    projector::Union{Nothing,P}
    M2C::M
    calibration::C
end

function ao_calibration(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS;
    n_modes::Int=size(dm.state.modes, 2), amplitude::Real=1e-9,
    projector::Bool=true, basis::Union{Nothing,ModalBasis}=nothing,
    method::KLBasisMethod=KLDMModes(), atm::Union{Nothing,AbstractAtmosphere}=nothing,
    build_backend::BuildBackend=default_build_backend(dm.state.coefs))

    if basis === nothing
        basis = modal_basis(dm, tel; n_modes=n_modes, projector=projector, method=method, atm=atm)
    end

    imat = interaction_matrix(dm, wfs, tel, basis.M2C; amplitude=amplitude)
    calib = CalibrationVault(imat.matrix; build_backend=build_backend)
    return AOCalibration(basis.basis, basis.projector, basis.M2C, calib)
end

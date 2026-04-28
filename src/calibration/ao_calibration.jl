#
# AO calibration assembly
#
# `ao_calibration` packages the standard modal calibration flow:
# 1. choose or build a modal basis for the DM
# 2. measure the WFS interaction matrix in that command basis
# 3. build a `ControlMatrix` inverse operator from the measured matrix
#
# The result is the compact object used by modal reconstructors and runtime
# controller setup.
#
struct AOCalibration{T<:AbstractFloat,
    B<:AbstractMatrix{T},
    P<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    C<:ControlMatrix{T}}
    basis::B
    projector::Union{Nothing,P}
    M2C::M
    calibration::C
end

@inline modal_to_command(calib::AOCalibration) = calib.M2C
@inline sampled_basis(calib::AOCalibration) = calib.basis
@inline modal_projector(calib::AOCalibration) = calib.projector
@inline control_matrix(calib::AOCalibration) = calib.calibration

"""
    ao_calibration(tel, dm, wfs; ...)

Build the standard modal AO calibration package.

If no basis is provided, the function first constructs one from the DM and
telescope, then measures the interaction matrix in that modal command basis and
builds a `ControlMatrix` inverse from the result.
"""
function ao_calibration(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS;
    n_modes::Int=size(dm.state.modes, 2), amplitude::Real=1e-9,
    projector::Bool=true, basis::Union{Nothing,ModalBasis}=nothing,
    method::KLBasisMethod=KLDMModes(), atm::Union{Nothing,AbstractAtmosphere}=nothing,
    build_backend::BuildBackend=default_runtime_calibration_build_backend(dm.state.coefs))

    if basis === nothing
        basis = modal_basis(dm, tel; n_modes=n_modes, projector=projector, method=method, atm=atm)
    end

    imat = interaction_matrix(dm, wfs, tel, basis.M2C; amplitude=amplitude)
    calib = ControlMatrix(imat.matrix; build_backend=build_backend)
    return AOCalibration(basis.basis, basis.projector, basis.M2C, calib)
end

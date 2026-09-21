#
# Modal basis construction
#
# This file builds the command-basis operators used for modal calibration and
# control.
#
# Supported constructions:
# - `InfluenceFunctionEigenbasis`: delegated sampled influence-function Gram
#   eigenbasis from AdaptiveOpticsCalibration
# - `KarhunenLoeveBasis`: finite-DM atmospheric covariance eigenbasis delegated
#   to AdaptiveOpticsCalibration after AOS constructs the physical covariance
#
# The resulting `M2C` operator maps modal coefficients to actuator commands,
# while `basis` stores the corresponding OPD modes on the pupil grid.
#
const _AOC_MODAL_BASES = AdaptiveOpticsCalibration.ModalBases

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

function ModalBasis(
    modal_to_command::M,
    basis::B,
    ::Nothing,
) where {T<:AbstractFloat,M<:AbstractMatrix{T},B<:AbstractMatrix{T}}
    return ModalBasis{T,M,B,M}(modal_to_command, basis, nothing)
end

@inline modal_to_command(basis::ModalBasis) = basis.M2C
@inline sampled_basis(basis::ModalBasis) = basis.basis
@inline modal_projector(basis::ModalBasis) = basis.projector

function dm_basis(dm::DeformableMirror, tel::Telescope)
    n = tel.params.resolution
    n_modes = size(dm.state.modes, 2)
    return reshape(sampled_influence_matrix(dm), n, n, n_modes)
end

function basis_from_m2c(dm::DeformableMirror, tel::Telescope, M2C::AbstractMatrix)
    n = tel.params.resolution
    basis_mat = sampled_influence_matrix(dm) * M2C
    return reshape(basis_mat, n, n, size(M2C, 2))
end

"""
    basis_projector(basis; method=..., build_backend=...)

Construct a projector from sampled basis vectors back into modal coefficients.

The configured Calibration method owns the inverse construction and its
rank-deficiency behavior. The accepted product is materialized on the selected
runtime backend.
"""
function basis_projector(basis::AbstractMatrix{T};
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse=
        _default_svd_inverse_method(T),
    build_backend::BuildBackend=
        default_runtime_calibration_build_backend(basis),
) where {T<:AbstractFloat}
    product = _prepare_svd_reconstructor(basis, method)
    return materialize_runtime_build_result(
        build_backend,
        basis,
        _AOC_RECONSTRUCTORS.reconstructor(product),
    )
end

function _modal_basis_projector(
    ::_AOC_MODAL_BASES.AbstractModalBasisMethod,
    basis::AbstractMatrix,
    ::Telescope,
)
    return basis_projector(basis)
end

function _modal_basis_projector(
    ::KarhunenLoeveBasis,
    basis::AbstractMatrix{T},
    tel::Telescope,
) where {T<:AbstractFloat}
    host_basis = prepare_build_matrix(CPUBuildBackend(), basis)
    support = vec(Array(pupil_mask(tel)))
    support_count = count(support)
    support_count > 0 || throw(InvalidConfiguration("pupil support must be nonempty"))
    weight = inv(T(support_count))
    host_projector = Matrix{T}(undef, size(host_basis, 2), size(host_basis, 1))
    @inbounds for sample in axes(host_basis, 1), mode in axes(host_basis, 2)
        host_projector[mode, sample] = support[sample] ?
            host_basis[sample, mode] * weight : zero(T)
    end
    return materialize_runtime_build_result(
        default_runtime_calibration_build_backend(basis),
        basis,
        host_projector,
    )
end

@inline _default_modal_mode_count(
    ::_AOC_MODAL_BASES.AbstractModalBasisMethod,
    dm::DeformableMirror,
) = size(dm.state.modes, 2)

@inline function _default_modal_mode_count(
    method::KarhunenLoeveBasis,
    dm::DeformableMirror,
)
    corrector_count = size(dm.state.modes, 2)
    return method.remove_piston ? max(corrector_count - 1, 1) : corrector_count
end

@inline _resolve_modal_mode_count(
    method::_AOC_MODAL_BASES.AbstractModalBasisMethod,
    dm::DeformableMirror,
    ::Nothing,
) = _default_modal_mode_count(method, dm)

@inline _resolve_modal_mode_count(
    ::_AOC_MODAL_BASES.AbstractModalBasisMethod,
    ::DeformableMirror,
    n_modes::Int,
) = n_modes

"""
    projected_atmospheric_opd_covariance(sampled_influences, pupil_support,
                                         resolution, sampling_m, atmosphere)

Construct the corrector-coordinate atmospheric OPD covariance used by the
finite-DM Karhunen–Loève calibration. The normalized pupil-supported influence
functions are transformed on a doubled periodic spectral grid. The von Kármán
phase spectrum uses cycles per metre and is converted to OPD at the
atmosphere's reference wavelength before frequency-cell quadrature.
"""
function projected_atmospheric_opd_covariance(
    sampled_influences::AbstractMatrix{T},
    pupil_support::AbstractVector{Bool},
    resolution::Int,
    sampling_m::T,
    atm::AbstractAtmosphere,
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(sampled_influences, pupil_support)
    sample_count = Base.Checked.checked_mul(resolution, resolution)
    size(sampled_influences, 1) == sample_count || throw(DimensionMismatchError(
        "sampled influences must have one row per pupil-grid sample",
    ))
    length(pupil_support) == sample_count || throw(DimensionMismatchError(
        "pupil support must have one entry per pupil-grid sample",
    ))
    support_count = count(pupil_support)
    support_count > 0 || throw(InvalidConfiguration("pupil support must be nonempty"))
    isfinite(sampling_m) && sampling_m > zero(T) || throw(InvalidConfiguration(
        "pupil sampling must be positive and finite",
    ))
    all(isfinite, sampled_influences) || throw(InvalidConfiguration(
        "sampled influences must contain only finite values",
    ))

    corrector_count = size(sampled_influences, 2)
    spectral_resolution = Base.Checked.checked_mul(2, resolution)
    freqs = Vector{T}(undef, spectral_resolution)
    fftfreq!(freqs, spectral_resolution; d=sampling_m)
    radial_frequency = Matrix{T}(undef, spectral_resolution, spectral_resolution)
    @inbounds for column in 1:spectral_resolution, row in 1:spectral_resolution
        radial_frequency[row, column] = hypot(freqs[row], freqs[column])
    end
    r0, L0, reference_wavelength_m = atmospheric_spectrum_parameters(atm)
    opd_scale = (T(reference_wavelength_m) / T(2π))^2
    opd_psd = phase_spectrum(radial_frequency, T(r0), T(L0))
    opd_psd .*= opd_scale

    transformed_influences = Matrix{Complex{T}}(
        undef,
        spectral_resolution * spectral_resolution,
        corrector_count,
    )
    buffer = zeros(Complex{T}, spectral_resolution, spectral_resolution)
    fft_plan = plan_fft_backend!(buffer)
    support_weight = inv(T(support_count))
    offset = (spectral_resolution - resolution) ÷ 2
    @inbounds for corrector in 1:corrector_count
        fill!(buffer, zero(Complex{T}))
        for column in 1:resolution, row in 1:resolution
            sample = (column - 1) * resolution + row
            if pupil_support[sample]
                buffer[offset + row, offset + column] = complex(
                    support_weight * sampled_influences[sample, corrector],
                    zero(T),
                )
            end
        end
        execute_fft_plan!(buffer, fft_plan)
        @views transformed_influences[:, corrector] .= reshape(buffer, :)
    end

    psd_vector = reshape(opd_psd, :)
    weighted = similar(transformed_influences)
    @inbounds for corrector in 1:corrector_count
        destination = view(weighted, :, corrector)
        transformed = view(transformed_influences, :, corrector)
        for index in eachindex(psd_vector)
            destination[index] = transformed[index] * psd_vector[index]
        end
    end
    frequency_step = inv(T(spectral_resolution) * sampling_m)
    covariance = real(adjoint(transformed_influences) * weighted)
    covariance .*= frequency_step^2
    @inbounds for column in axes(covariance, 2)
        for row in 1:column
            value = (covariance[row, column] + covariance[column, row]) / T(2)
            covariance[row, column] = value
            covariance[column, row] = value
        end
    end
    return covariance
end

function atmospheric_spectrum_parameters(atm::KolmogorovAtmosphere)
    return atm.params.r0, atm.params.L0, atm.params.reference_wavelength_m
end

function atmospheric_spectrum_parameters(atm::MultiLayerAtmosphere)
    return atm.params.r0, atm.params.L0, atm.params.reference_wavelength_m
end

function atmospheric_spectrum_parameters(atm::AbstractAtmosphere)
    throw(InvalidConfiguration(
        "atmospheric spectrum parameters are not defined for $(typeof(atm))",
    ))
end

"""
    modal_basis(dm, tel; ...)

Build the modal command basis used by AO calibration and control.

This returns both the modal-to-command operator and the sampled pupil-space
basis, with an optional projector back into modal coordinates.
"""
function modal_basis(dm::DeformableMirror, tel::Telescope;
    n_modes::Union{Nothing,Int}=nothing,
    projector::Bool=true,
    method::_AOC_MODAL_BASES.AbstractModalBasisMethod=
        _AOC_MODAL_BASES.InfluenceFunctionEigenbasis(),
    atm::Union{Nothing,AbstractAtmosphere}=nothing)
    resolved_mode_count = _resolve_modal_mode_count(method, dm, n_modes)
    M2C, basis = modal_basis_components(
        method,
        dm,
        tel,
        atm;
        n_modes=resolved_mode_count,
    )
    basis_mat = reshape(basis, :, size(basis, 3))
    proj = projector ? _modal_basis_projector(method, basis_mat, tel) : nothing
    return ModalBasis(M2C, basis_mat, proj)
end

function _materialize_modal_basis(product, sampled_influences::AbstractMatrix,
    tel::Telescope, n_modes::Int)
    build_backend = default_runtime_calibration_build_backend(sampled_influences)
    modal_to_command = materialize_runtime_build_result(
        build_backend,
        sampled_influences,
        _AOC_MODAL_BASES.modal_to_command(product),
    )
    sampled_modes = materialize_runtime_build_result(
        build_backend,
        sampled_influences,
        _AOC_MODAL_BASES.sampled_modes(product),
    )
    basis = reshape(
        sampled_modes,
        tel.params.resolution,
        tel.params.resolution,
        n_modes,
    )
    return modal_to_command, basis
end

function modal_basis_components(
    method::_AOC_MODAL_BASES.InfluenceFunctionEigenbasis,
    dm::DeformableMirror,
    tel::Telescope,
    ::Union{Nothing,AbstractAtmosphere};
    n_modes::Int,
)
    sampled_influences = sampled_influence_matrix(dm)
    T = eltype(sampled_influences)
    host_influences = prepare_build_matrix(CPUBuildBackend(), sampled_influences)
    support = vec(Array(pupil_mask(tel)))
    specification = _AOC_MODAL_BASES.SampledInfluenceBasisSpecification(
        size(host_influences, 1),
        size(host_influences, 2),
        n_modes,
        support,
        T,
    )
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    product = AdaptiveOpticsCalibration.process(plan, host_influences)
    return _materialize_modal_basis(product, sampled_influences, tel, n_modes)
end

function modal_basis_components(method::KarhunenLoeveBasis, dm::DeformableMirror,
    tel::Telescope, ::Nothing;
    n_modes::Int)
    throw(InvalidConfiguration("KarhunenLoeveBasis requires an atmosphere"))
end

function modal_basis_components(method::KarhunenLoeveBasis, dm::DeformableMirror,
    tel::Telescope, atm::AbstractAtmosphere;
    n_modes::Int)
    sampled_influences = sampled_influence_matrix(dm)
    host_influences = prepare_build_matrix(CPUBuildBackend(), sampled_influences)
    support = vec(Array(pupil_mask(tel)))
    sampling_x, sampling_y = tel.aperture.sampling_m
    sampling_x == sampling_y || throw(InvalidConfiguration(
        "Karhunen–Loève spectral covariance requires square pupil sampling",
    ))
    projected_covariance = projected_atmospheric_opd_covariance(
        host_influences,
        support,
        tel.params.resolution,
        eltype(host_influences)(sampling_x),
        atm,
    )
    specification = _AOC_MODAL_BASES.KarhunenLoeveBasisSpecification(
        size(host_influences, 1),
        size(host_influences, 2),
        n_modes,
        support,
        eltype(host_influences),
    )
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    product = AdaptiveOpticsCalibration.process(
        plan,
        _AOC_MODAL_BASES.KarhunenLoeveBasisInputs(
            host_influences,
            projected_covariance,
        ),
    )
    return _materialize_modal_basis(product, sampled_influences, tel, n_modes)
end

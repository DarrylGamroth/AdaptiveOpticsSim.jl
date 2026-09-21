# Cold bridge from a plant response matrix to the reusable Calibration product.
const _AOC_RECONSTRUCTORS = AdaptiveOpticsCalibration.Reconstructors

@inline _default_svd_inverse_method(::Type{T}) where {T<:AbstractFloat} =
    _AOC_RECONSTRUCTORS.TSVDInverse(rtol=sqrt(eps(T)))

function _prepare_svd_reconstructor(
    interaction_matrix::AbstractMatrix{T},
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse,
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(interaction_matrix)
    specification = _AOC_RECONSTRUCTORS.ReconstructorSpecification(
        size(interaction_matrix, 1),
        size(interaction_matrix, 2),
        T,
    )
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    inputs = _AOC_RECONSTRUCTORS.SVDReconstructorInputs(
        Matrix{T}(interaction_matrix),
    )
    return AdaptiveOpticsCalibration.process(plan, inputs)
end

# Prepared modal OPD expansion

struct _OwnedModalOPDExpansionPlan end
const _OWNED_MODAL_OPD_EXPANSION_PLAN = _OwnedModalOPDExpansionPlan()

"""
    ModalOPDExpansionPlan(basis, pupil_support)

Prepare a run-immutable modal OPD expansion contract.

`basis[:, :, k]` is a dimensionless pupil-plane OPD mode. The corresponding
coefficient is expressed in metres. Construction snapshots both arrays on
their existing compute device. The plan owns no persistent state or scratch
workspace.
"""
struct ModalOPDExpansionPlan{
    T<:AbstractFloat,
    B<:AbstractArray{T,3},
    P<:AbstractMatrix{Bool},
}
    basis::B
    pupil_support::P

    function ModalOPDExpansionPlan(
        ::_OwnedModalOPDExpansionPlan,
        basis::B,
        pupil_support::P,
    ) where {T<:AbstractFloat,B<:AbstractArray{T,3},P<:AbstractMatrix{Bool}}
        Base.require_one_based_indexing(basis, pupil_support)
        axes(basis, 1) == axes(pupil_support, 1) &&
            axes(basis, 2) == axes(pupil_support, 2) || throw(
            DimensionMismatchError(
                "modal OPD basis pupil axes must match pupil-support axes",
            ),
        )
        compute_device(basis) == compute_device(pupil_support) || throw(
            InvalidConfiguration(
                "modal OPD basis and pupil support must occupy one compute device",
            ),
        )
        return new{T,B,P}(basis, pupil_support)
    end
end

function ModalOPDExpansionPlan(
    basis::AbstractArray{T,3},
    pupil_support::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(basis, pupil_support)
    return ModalOPDExpansionPlan(
        _OWNED_MODAL_OPD_EXPANSION_PLAN,
        copy(basis),
        copy(pupil_support),
    )
end

function _replace_modal_opd_basis(
    plan::ModalOPDExpansionPlan{T},
    basis::AbstractArray{T,3},
) where {T}
    axes(basis) == axes(plan.basis) || throw(DimensionMismatchError(
        "replacement modal OPD basis axes must match the prepared basis axes",
    ))
    return ModalOPDExpansionPlan(
        _OWNED_MODAL_OPD_EXPANSION_PLAN,
        copy(basis),
        plan.pupil_support,
    )
end

function _replace_modal_opd_pupil_support(
    plan::ModalOPDExpansionPlan,
    pupil_support::AbstractMatrix{Bool},
)
    axes(pupil_support) == axes(plan.pupil_support) || throw(
        DimensionMismatchError(
            "replacement pupil-support axes must match the prepared axes",
        ),
    )
    return ModalOPDExpansionPlan(
        _OWNED_MODAL_OPD_EXPANSION_PLAN,
        plan.basis,
        copy(pupil_support),
    )
end

"""
    combine_basis!(opd, plan::ModalOPDExpansionPlan, coefficients)

Evaluate a prepared modal OPD expansion into caller-owned `opd` storage.

The output and coefficients must use one-based axes, match the plan exactly,
and occupy the plan's compute device. Successful repeated calls allocate no
Julia heap storage after warmup for ordinary CPU arrays.
"""
function combine_basis!(
    opd::AbstractMatrix{T},
    plan::ModalOPDExpansionPlan{T},
    coefficients::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(opd, coefficients)
    axes(opd, 1) == axes(plan.pupil_support, 1) &&
        axes(opd, 2) == axes(plan.pupil_support, 2) || throw(
        DimensionMismatchError(
            "modal OPD output axes must match the prepared pupil-support axes",
        ),
    )
    axes(coefficients, 1) == axes(plan.basis, 3) || throw(
        DimensionMismatchError(
            "modal coefficients must match the prepared basis-mode axis",
        ),
    )
    device = compute_device(opd)
    device == compute_device(coefficients) &&
        device == compute_device(plan.basis) &&
        device == compute_device(plan.pupil_support) || throw(
        InvalidConfiguration(
            "modal OPD input, output, basis, and pupil support must occupy one compute device",
        ),
    )
    return combine_basis!(opd, plan.basis, coefficients, plan.pupil_support)
end

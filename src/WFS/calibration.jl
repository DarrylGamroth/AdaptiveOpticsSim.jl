@inline calibration_wavelength(src::AbstractSource, ::Type{T}) where {T<:AbstractFloat} = T(wavelength(src))

@inline calibration_signature(src::AbstractSource) = source_measurement_signature(src)

@inline calibration_matches(calibrated::Bool, stored_λ, λ) = calibrated && stored_λ == λ

@inline calibration_matches(calibrated::Bool, stored_λ, λ, stored_sig::UInt, sig::UInt) =
    calibrated && stored_λ == λ && stored_sig == sig

function save_zero_opd!(tel::Telescope)
    saved = copy(tel.state.opd)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    return saved
end

@inline function restore_opd!(tel::Telescope, saved_opd::AbstractMatrix)
    copyto!(tel.state.opd, saved_opd)
    return tel
end

@inline function store_reference_signal!(reference::AbstractMatrix, signal::AbstractMatrix, slopes::AbstractVector)
    copyto!(reference, signal)
    fill!(slopes, zero(eltype(slopes)))
    return reference
end

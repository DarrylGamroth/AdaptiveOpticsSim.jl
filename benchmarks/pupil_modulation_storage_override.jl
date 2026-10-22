# Diagnostic preparation-only storage alternatives. Never included by package code.
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
import CUDA
const PupilStorageOptics = AdaptiveOpticsSim.Optics
const PupilStorageWFS = AdaptiveOpticsSim.WavefrontSensors
const PUPIL_STORAGE = get(ENV,"AOS_PUPIL_STORAGE","separable")
PUPIL_STORAGE in ("dense","separable","computed") || error("invalid pupil storage diagnostic")

struct ComputedPupilModulation{P,C,O,V}
    policy::P
    coordinates::C
    offsets::O
    amplitude_weights::V
end
PupilStorageOptics.modulation_point_count(m::ComputedPupilModulation) = size(m.offsets,2)
PupilStorageWFS._pyramid_modulation_resolution(m::ComputedPupilModulation) = length(m.coordinates)
PupilStorageOptics._modulation_phase_storage(m::ComputedPupilModulation) = (m.coordinates,m.offsets)
@inline function PupilStorageOptics._modulation_phase(storage::Tuple{C,O},i,j,p) where {
    C<:AbstractVector{<:Real}, O<:AbstractMatrix{<:Real}}
    @inbounds return cis(storage[2][1,p]*storage[1][j]+storage[2][2,p]*storage[1][i])
end
function PupilStorageWFS._require_exact_focal_plane_modulation_target(
    m::ComputedPupilModulation,target::Backends.AbstractComputeDevice,label::AbstractString)
    PupilStorageWFS._require_exact_wfs_storage_target(m.coordinates,target,:wfs_optics,"$label coordinates")
    PupilStorageWFS._require_exact_wfs_storage_target(m.offsets,target,:wfs_optics,"$label offsets")
    return m
end
function PupilStorageWFS._prepare_pyramid_modulation(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend},policy,resolution,field,::Type{T}) where {T}
    PUPIL_STORAGE == "dense" && return PupilStorageOptics.prepare_focal_plane_modulation(policy,resolution,field,T)
    PUPIL_STORAGE == "separable" && return PupilStorageOptics.prepare_separable_focal_plane_modulation(policy,resolution,field,T)
    count=PupilStorageOptics.modulation_point_count(policy)
    coordinates=similar(field,T,resolution)
    offsets=similar(field,T,2,count)
    copyto!(coordinates,collect(range(-T(pi),T(pi);length=resolution)))
    host_offsets=T[PupilStorageOptics.modulation_offset(policy,p,T)[axis] for axis in 1:2,p in 1:count]
    copyto!(offsets,host_offsets)
    weights=PupilStorageOptics.FixedSizeVectorDefault{T}([
        sqrt(PupilStorageOptics.modulation_weight(policy,p,T)) for p in 1:count])
    return ComputedPupilModulation(policy,coordinates,offsets,weights)
end

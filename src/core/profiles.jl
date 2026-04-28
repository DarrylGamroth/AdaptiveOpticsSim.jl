abstract type FidelityProfile end

struct ScientificProfile <: FidelityProfile end
struct FastProfile <: FidelityProfile end

struct ProfileBundle{
    A<:FidelityProfile,
    C<:FidelityProfile,
    D<:FidelityProfile,
    L<:FidelityProfile,
    T<:FidelityProfile,
} <: FidelityProfile
    atmosphere::A
    calibration::C
    detector::D
    lift::L
    tomography::T
end

default_fidelity_profile() = ScientificProfile()

function ProfileBundle(base::FidelityProfile=default_fidelity_profile();
    atmosphere::FidelityProfile=base,
    calibration::FidelityProfile=base,
    detector::FidelityProfile=base,
    lift::FidelityProfile=base,
    tomography::FidelityProfile=base)
    return ProfileBundle(atmosphere, calibration, detector, lift, tomography)
end

atmosphere_profile(profile::FidelityProfile) = profile
calibration_profile(profile::FidelityProfile) = profile
detector_profile(profile::FidelityProfile) = profile
lift_profile(profile::FidelityProfile) = profile
tomography_profile(profile::FidelityProfile) = profile

atmosphere_profile(profile::ProfileBundle) = profile.atmosphere
calibration_profile(profile::ProfileBundle) = profile.calibration
detector_profile(profile::ProfileBundle) = profile.detector
lift_profile(profile::ProfileBundle) = profile.lift
tomography_profile(profile::ProfileBundle) = profile.tomography

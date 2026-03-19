abstract type FidelityProfile end

struct ScientificProfile <: FidelityProfile end
struct FastProfile <: FidelityProfile end

default_fidelity_profile() = ScientificProfile()


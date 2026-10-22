# Benchmark-only preparation override. Include before constructing any sensor.
# Zero selects the unchanged baseline; -1 selects the production policy.
# No production environment knob is introduced.
using AdaptiveOpticsSim
const ModulationBatchWFS = AdaptiveOpticsSim.WavefrontSensors
const MODULATION_BATCH = Ref(0)
@eval ModulationBatchWFS function _pyramid_modulation_batch_size(points::Int, bytes::Int,
    limit::Int=_PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT)
    points >= 1 && bytes >= 1 || throw(InvalidConfiguration("invalid modulation batch dimensions"))
    requested = Main.MODULATION_BATCH[]
    bound = max(1, div(_PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES, bytes))
    if requested > 0
        1 <= requested <= bound && points % requested == 0 ||
            throw(InvalidConfiguration("benchmark batch must divide points and fit scratch budget"))
        return requested
    end
    requested >= -1 || throw(InvalidConfiguration("invalid benchmark batch override"))
    selected_limit = requested == -1 ? limit : _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT
    for candidate in min(points, selected_limit, bound):-1:1
        points % candidate == 0 && return candidate
    end
    return 1
end

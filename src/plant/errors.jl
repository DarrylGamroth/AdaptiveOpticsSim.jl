"""Invalid canonical plant timestamp, duration, conversion, or arithmetic."""
struct PlantTimeError <: AdaptiveOpticsSimError
    operation::Symbol
    reason::Symbol
    msg::String
end

"""Invalid periodic schedule, event generator, or deterministic scheduler operation."""
struct PlantScheduleError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Invalid immutable plant, optical-path, or acquisition declaration."""
struct PlantDefinitionError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Invalid plant-command schema, payload, admission, or lifecycle transition."""
struct PlantCommandError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

"""Invalid preparation or prepared binding for a plant component."""
struct PlantPreparationError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Invalid prepared detector-acquisition event or lifecycle transition."""
struct DetectorAcquisitionError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

struct RuntimeOutputRequirements
    slopes::Bool
    wfs_pixels::Bool
    science_pixels::Bool
end

RuntimeOutputRequirements(; slopes::Bool=true, wfs_pixels::Bool=false, science_pixels::Bool=false) =
    RuntimeOutputRequirements(slopes, wfs_pixels, science_pixels)

struct GroupedRuntimeOutputRequirements
    wfs_frames::Bool
    science_frames::Bool
    wfs_stack::Bool
    science_stack::Bool
end

GroupedRuntimeOutputRequirements(; wfs_frames::Bool=true, science_frames::Bool=true,
    wfs_stack::Bool=true, science_stack::Bool=true) =
    GroupedRuntimeOutputRequirements(wfs_frames, science_frames, wfs_stack, science_stack)

struct RuntimeOutputPlan
    slopes::Bool
    wfs_frame::Bool
    science_frame::Bool
end

struct GroupedRuntimeOutputPlan
    wfs_frames::Bool
    science_frames::Bool
    wfs_stack::Bool
    science_stack::Bool
end

default_runtime_outputs(; wfs_detector=nothing, science_detector=nothing) =
    RuntimeOutputRequirements(slopes=true, wfs_pixels=!isnothing(wfs_detector), science_pixels=!isnothing(science_detector))

default_grouped_runtime_outputs(interfaces...) =
    GroupedRuntimeOutputRequirements(
        wfs_frames=any(i -> !isnothing(simulation_wfs_frame(i)), interfaces),
        science_frames=any(i -> !isnothing(simulation_science_frame(i)), interfaces),
        wfs_stack=true,
        science_stack=true,
    )

@inline runtime_output_plan(outputs::RuntimeOutputRequirements, wfs_detector, science_detector) =
    RuntimeOutputPlan(outputs.slopes, outputs.wfs_pixels && !isnothing(wfs_detector), outputs.science_pixels && !isnothing(science_detector))

@inline runtime_output_plan(runtime) = runtime_output_plan(runtime.outputs, runtime.wfs_detector, runtime.science_detector)

@inline function grouped_runtime_output_plan(outputs::GroupedRuntimeOutputRequirements, interfaces::Tuple)
    wfs_frames = outputs.wfs_frames && any(i -> !isnothing(simulation_wfs_frame(i)), interfaces)
    science_frames = outputs.science_frames && any(i -> !isnothing(simulation_science_frame(i)), interfaces)
    wfs_stack = outputs.wfs_stack && _supports_grouped_frame_stack(map(simulation_wfs_frame, interfaces))
    science_stack = outputs.science_stack && _supports_grouped_frame_stack(map(simulation_science_frame, interfaces))
    return GroupedRuntimeOutputPlan(wfs_frames, science_frames, wfs_stack, science_stack)
end

@inline grouped_runtime_output_plan(interface) =
    grouped_runtime_output_plan(interface.outputs, interface.interfaces)

@inline requires_runtime_slopes(runtime) = runtime_output_plan(runtime).slopes
@inline requires_runtime_wfs_pixels(runtime) = runtime_output_plan(runtime).wfs_frame
@inline requires_runtime_science_pixels(runtime) = runtime_output_plan(runtime).science_frame

@inline runtime_outputs(runtime) = runtime.outputs
@inline grouped_runtime_outputs(interface) = interface.outputs

@inline requires_grouped_wfs_frames(interface) = grouped_runtime_output_plan(interface).wfs_frames
@inline requires_grouped_science_frames(interface) = grouped_runtime_output_plan(interface).science_frames
@inline requires_grouped_wfs_stack(interface) = grouped_runtime_output_plan(interface).wfs_stack
@inline requires_grouped_science_stack(interface) = grouped_runtime_output_plan(interface).science_stack

@inline function _grouped_frame_signature(frame)
    return (size(frame), eltype(frame), ndims(frame))
end

function _supports_grouped_frame_stack(frames)
    reference = nothing
    @inbounds for frame in frames
        isnothing(frame) && return false
        sig = _grouped_frame_signature(frame)
        if isnothing(reference)
            reference = sig
        elseif sig != reference
            return false
        end
    end
    return !isnothing(reference)
end

@inline _grouped_stack_dims(frame, count::Int) = (size(frame)..., count)

function grouped_frame_stack_buffer(frames)
    _supports_grouped_frame_stack(frames) || return nothing
    reference = first(frames)
    return similar(reference, _grouped_stack_dims(reference, length(frames)))
end

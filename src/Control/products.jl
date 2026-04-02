struct RuntimeProductRequirements
    slopes::Bool
    wfs_pixels::Bool
    science_pixels::Bool
end

RuntimeProductRequirements(; slopes::Bool=true, wfs_pixels::Bool=false, science_pixels::Bool=false) =
    RuntimeProductRequirements(slopes, wfs_pixels, science_pixels)

struct GroupedRuntimeProductRequirements
    wfs_frames::Bool
    science_frames::Bool
    wfs_stack::Bool
    science_stack::Bool
end

GroupedRuntimeProductRequirements(; wfs_frames::Bool=true, science_frames::Bool=true,
    wfs_stack::Bool=true, science_stack::Bool=true) =
    GroupedRuntimeProductRequirements(wfs_frames, science_frames, wfs_stack, science_stack)

struct RuntimeProductPlan
    slopes::Bool
    wfs_frame::Bool
    science_frame::Bool
end

struct GroupedRuntimeProductPlan
    wfs_frames::Bool
    science_frames::Bool
    wfs_stack::Bool
    science_stack::Bool
end

default_runtime_products(; wfs_detector=nothing, science_detector=nothing) =
    RuntimeProductRequirements(slopes=true, wfs_pixels=!isnothing(wfs_detector), science_pixels=!isnothing(science_detector))

default_grouped_runtime_products(interfaces...) =
    GroupedRuntimeProductRequirements(
        wfs_frames=any(i -> !isnothing(simulation_wfs_frame(i)), interfaces),
        science_frames=any(i -> !isnothing(simulation_science_frame(i)), interfaces),
        wfs_stack=true,
        science_stack=true,
    )

@inline runtime_product_plan(products::RuntimeProductRequirements, wfs_detector, science_detector) =
    RuntimeProductPlan(products.slopes, products.wfs_pixels && !isnothing(wfs_detector), products.science_pixels && !isnothing(science_detector))

@inline runtime_product_plan(runtime) = runtime_product_plan(runtime.products, runtime.wfs_detector, runtime.science_detector)

@inline function grouped_runtime_product_plan(products::GroupedRuntimeProductRequirements, interfaces::Tuple)
    wfs_frames = products.wfs_frames && any(i -> !isnothing(simulation_wfs_frame(i)), interfaces)
    science_frames = products.science_frames && any(i -> !isnothing(simulation_science_frame(i)), interfaces)
    wfs_stack = products.wfs_stack && _supports_grouped_frame_stack(map(simulation_wfs_frame, interfaces))
    science_stack = products.science_stack && _supports_grouped_frame_stack(map(simulation_science_frame, interfaces))
    return GroupedRuntimeProductPlan(wfs_frames, science_frames, wfs_stack, science_stack)
end

@inline grouped_runtime_product_plan(interface) =
    grouped_runtime_product_plan(interface.products, interface.interfaces)

@inline requires_runtime_slopes(runtime) = runtime_product_plan(runtime).slopes
@inline requires_runtime_wfs_pixels(runtime) = runtime_product_plan(runtime).wfs_frame
@inline requires_runtime_science_pixels(runtime) = runtime_product_plan(runtime).science_frame

@inline runtime_products(runtime) = runtime.products
@inline grouped_runtime_products(interface) = interface.products

@inline requires_grouped_wfs_frames(interface) = grouped_runtime_product_plan(interface).wfs_frames
@inline requires_grouped_science_frames(interface) = grouped_runtime_product_plan(interface).science_frames
@inline requires_grouped_wfs_stack(interface) = grouped_runtime_product_plan(interface).wfs_stack
@inline requires_grouped_science_stack(interface) = grouped_runtime_product_plan(interface).science_stack

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

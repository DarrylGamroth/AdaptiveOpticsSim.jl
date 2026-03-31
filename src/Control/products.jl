struct RuntimeProductRequirements
    slopes::Bool
    wfs_pixels::Bool
    science_pixels::Bool
end

RuntimeProductRequirements(; slopes::Bool=true, wfs_pixels::Bool=false, science_pixels::Bool=false) =
    RuntimeProductRequirements(slopes, wfs_pixels, science_pixels)

struct RuntimeProductPlan
    slopes::Bool
    wfs_frame::Bool
    science_frame::Bool
end

default_runtime_products(; wfs_detector=nothing, science_detector=nothing) =
    RuntimeProductRequirements(slopes=true, wfs_pixels=!isnothing(wfs_detector), science_pixels=!isnothing(science_detector))

@inline runtime_product_plan(products::RuntimeProductRequirements, wfs_detector, science_detector) =
    RuntimeProductPlan(products.slopes, products.wfs_pixels && !isnothing(wfs_detector), products.science_pixels && !isnothing(science_detector))

@inline runtime_product_plan(runtime) = runtime_product_plan(runtime.products, runtime.wfs_detector, runtime.science_detector)

@inline requires_runtime_slopes(runtime) = runtime_product_plan(runtime).slopes
@inline requires_runtime_wfs_pixels(runtime) = runtime_product_plan(runtime).wfs_frame
@inline requires_runtime_science_pixels(runtime) = runtime_product_plan(runtime).science_frame

@inline runtime_products(runtime) = runtime.products

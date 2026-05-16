using AdaptiveOpticsSim
using BenchmarkTools

const MISREG_AD_FIELDS = (:influence_width, :actuator_x, :actuator_y, :shift_x, :shift_y,
    :rotation_deg, :radial_scaling, :tangential_scaling)

const MISREG_AD_EPSILONS = Dict(
    :influence_width => 1e-6,
    :actuator_x => 1e-6,
    :actuator_y => 1e-6,
    :shift_x => 1e-6,
    :shift_y => 1e-6,
    :rotation_deg => 1e-5,
    :radial_scaling => 1e-6,
    :tangential_scaling => 1e-6,
)

function finite_difference_gaussian_dm_influence_jacobian(tel, dm, actuator_index::Integer,
    field_order::Tuple; epsilons=MISREG_AD_EPSILONS)
    n = tel.params.resolution
    jacobian = Matrix{eltype(dm.state.modes)}(undef, n * n, length(field_order))
    @inbounds for (idx, field) in pairs(field_order)
        jacobian[:, idx] .= finite_difference_gaussian_dm_influence_column(
            tel, dm, actuator_index, field, epsilons[field])
    end
    return jacobian
end

function finite_difference_gaussian_dm_influence_column(tel, dm, actuator_index::Integer,
    field::Symbol, epsilon::Real)
    base = gaussian_dm_influence_field_value(dm, actuator_index, field)
    plus = gaussian_dm_influence_column_with_field(tel, dm, actuator_index, field, base + epsilon)
    minus = gaussian_dm_influence_column_with_field(tel, dm, actuator_index, field, base - epsilon)
    return (plus .- minus) ./ (2 * epsilon)
end

function gaussian_dm_influence_field_value(dm, actuator_index::Integer, field::Symbol)
    coords = AdaptiveOpticsSim.actuator_coordinates(AdaptiveOpticsSim.topology(dm))
    mis = dm.params.misregistration
    if field === :influence_width
        return influence_width(dm)
    elseif field === :actuator_x
        return coords[1, actuator_index]
    elseif field === :actuator_y
        return coords[2, actuator_index]
    elseif field === :shift_x
        return mis.shift_x
    elseif field === :shift_y
        return mis.shift_y
    elseif field === :rotation_deg
        return AdaptiveOpticsSim.rotation_deg(mis)
    elseif field === :radial_scaling
        return mis.radial_scaling
    elseif field === :tangential_scaling
        return mis.tangential_scaling
    end
    throw(ArgumentError("unsupported finite-difference field $(field)"))
end

function gaussian_dm_influence_column_with_field(tel, dm, actuator_index::Integer,
    field::Symbol, value::Real)
    T = eltype(dm.state.modes)
    mis = dm.params.misregistration
    topology = AdaptiveOpticsSim.topology(dm)
    model = influence_model(dm)
    if field === :influence_width
        perturbed = DeformableMirror(tel; topology=topology, influence_width=value,
            misregistration=mis, T=T)
    elseif field === :actuator_x || field === :actuator_y
        coords = Matrix(AdaptiveOpticsSim.actuator_coordinates(topology))
        coords[field === :actuator_x ? 1 : 2, actuator_index] = T(value)
        perturbed_topology = SampledActuatorTopology(coords; T=T)
        perturbed = DeformableMirror(tel; topology=perturbed_topology,
            influence_width=influence_width(dm), misregistration=mis, T=T)
    else
        perturbed_mis = AdaptiveOpticsSim.update_misregistration(mis, field, value)
        perturbed = DeformableMirror(tel; topology=topology, influence_model=model,
            misregistration=perturbed_mis, T=T)
    end
    return copy(@view perturbed.state.modes[:, actuator_index])
end

function bench_gaussian_dm_influence_ad()
    mis = Misregistration(shift_x=0.03, shift_y=-0.02, rotation_deg=1.5,
        radial_scaling=1.01, tangential_scaling=0.98)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.35, misregistration=mis)
    actuator_index = 6
    fields = MISREG_AD_FIELDS

    ad = AdaptiveOpticsSim.gaussian_dm_influence_parameter_jacobian(tel, dm, actuator_index;
        field_order=fields)
    fd = finite_difference_gaussian_dm_influence_jacobian(tel, dm, actuator_index, fields)
    max_abs_error = maximum(abs, ad.jacobian .- fd)
    @assert isapprox(ad.jacobian, fd; rtol=3e-4, atol=3e-7)

    AdaptiveOpticsSim.gaussian_dm_influence_parameter_jacobian(tel, dm, actuator_index;
        field_order=fields)
    finite_difference_gaussian_dm_influence_jacobian(tel, dm, actuator_index, fields)
    alloc_ad = @allocated AdaptiveOpticsSim.gaussian_dm_influence_parameter_jacobian(
        tel, dm, actuator_index; field_order=fields)
    alloc_fd = @allocated finite_difference_gaussian_dm_influence_jacobian(
        tel, dm, actuator_index, fields)

    println("Gaussian DM influence AD sensitivity")
    println("  fields: $(fields)")
    println("  max abs AD-vs-FD error: $(max_abs_error)")
    println("  AD allocations after warmup: $(alloc_ad) bytes")
    println("  FD allocations after warmup: $(alloc_fd) bytes")
    println("  AD benchmark:")
    display(@benchmark AdaptiveOpticsSim.gaussian_dm_influence_parameter_jacobian(
        $tel, $dm, $actuator_index; field_order=$fields))
    println("  finite-difference benchmark:")
    display(@benchmark finite_difference_gaussian_dm_influence_jacobian(
        $tel, $dm, $actuator_index, $fields))
end

bench_gaussian_dm_influence_ad()

function bench_misregistration_wfs_ad_probe()
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    basis = modal_basis(dm, tel; n_modes=2).M2C[:, 1:2]
    fields = [:shift_x, :shift_y, :rotation_deg]

    fd = AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields,
        sensitivity=:finite_difference)
    ad = AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields)
    max_abs_error = maximum(abs, ad.meta.D .- fd.meta.D)
    @assert isapprox(ad.meta.D, fd.meta.D; rtol=2e-3, atol=1e-9)

    AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields,
        sensitivity=:finite_difference)
    AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields)
    alloc_fd = @allocated AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields,
        sensitivity=:finite_difference)
    alloc_ad = @allocated AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis; n_mis_reg=length(fields), field_order=fields)

    println("Misregistration WFS-path default AD")
    println("  fields: $(Tuple(fields))")
    println("  max abs AD-vs-FD error: $(max_abs_error)")
    println("  AD allocations after warmup: $(alloc_ad) bytes")
    println("  FD allocations after warmup: $(alloc_fd) bytes")
    println("  AD benchmark:")
    display(@benchmark AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        $tel, $dm, $wfs, $basis; n_mis_reg=$(length(fields)), field_order=$fields))
    println("  finite-difference benchmark:")
    display(@benchmark AdaptiveOpticsSim.compute_meta_sensitivity_matrix(
        $tel, $dm, $wfs, $basis; n_mis_reg=$(length(fields)), field_order=$fields,
        sensitivity=:finite_difference))
end

bench_misregistration_wfs_ad_probe()

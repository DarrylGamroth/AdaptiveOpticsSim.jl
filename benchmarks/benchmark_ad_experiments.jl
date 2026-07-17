using AdaptiveOpticsSim
using BenchmarkTools

function bench_misregistration_wfs_ad_probe()
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
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

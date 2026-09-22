include(joinpath(@__DIR__, "common.jl"))

"""Form a geometric-wavefront reference response; this is not an RTC calibration."""
function geometric_wavefront_response_matrix(tel::Telescope, dm::DeformableMirror,
    sensor::ShackHartmannWFS, commands::AbstractMatrix{<:Real})
    layout = subaperture_layout(sensor.front_end)
    n_signal = 2 * length(layout.valid_mask)
    response = zeros(eltype(dm.state.coefs), n_signal, size(commands, 2))
    pupil = PupilFunction(tel)
    command = similar(dm.state.coefs)
    sampling_m = (tel.params.diameter / tel.params.resolution,
        tel.params.diameter / tel.params.resolution)

    @inbounds for mode in axes(commands, 2)
        copyto!(command, view(commands, :, mode))
        set_command!(dm, command)
        update_surface!(dm)
        reset_opd!(pupil)
        apply_surface!(pupil, dm, DMReplace())
        geometric_wavefront_slopes!(view(response, :, mode), pupil.opd,
            layout.valid_mask, sampling_m)
    end
    return response
end

function main(; resolution::Int=16)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=3, influence_width=0.35)
    sensor = ShackHartmannWFS(tel; n_lenslets=2)
    commands = modal_basis(dm, tel; n_modes=3).M2C
    fields = (:shift_x, :shift_y)
    epsilon = Misregistration(shift_x=1e-4, shift_y=1e-4, T=Float64)

    reference = geometric_wavefront_response_matrix(tel, dm, sensor, commands)
    sensitivity_columns = Matrix{Float64}(undef, length(reference),
        length(fields))
    for (index, field) in enumerate(fields)
        δ = Calibration.misregistration_component(epsilon, field)
        dm_plus = DeformableMirror(tel; n_act=dm.params.n_act,
            influence_model=influence_model(dm),
            misregistration=Calibration.update_misregistration(
                Misregistration(T=Float64), field, δ))
        dm_minus = DeformableMirror(tel; n_act=dm.params.n_act,
            influence_model=influence_model(dm),
            misregistration=Calibration.update_misregistration(
                Misregistration(T=Float64), field, -δ))
        plus = geometric_wavefront_response_matrix(tel, dm_plus, sensor, commands)
        minus = geometric_wavefront_response_matrix(tel, dm_minus, sensor, commands)
        sensitivity_columns[:, index] .= vec((plus .- minus) ./ (2δ))
    end

    meta = Calibration.MetaSensitivity(
        Calibration.ControlMatrix(sensitivity_columns),
        Calibration.ControlMatrix(reference; invert=false),
        epsilon,
        collect(fields),
    )
    sprint = Calibration.SPRINT(meta, Misregistration(T=Float64),
        Misregistration(T=Float64), false, :finite_difference)
    injected = Misregistration(shift_x=5e-4, shift_y=-5e-4, T=Float64)
    dm_in = DeformableMirror(tel; n_act=dm.params.n_act,
        influence_model=influence_model(dm), misregistration=injected)
    response_in = geometric_wavefront_response_matrix(tel, dm_in, sensor, commands)
    estimate = Calibration.estimate!(sprint, response_in; precision=4)

    @info "SPRINT geometric-reference tutorial complete" shift_x=estimate.shift_x shift_y=estimate.shift_y
    return (
        injected=injected,
        estimate=estimate,
        response=reference,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

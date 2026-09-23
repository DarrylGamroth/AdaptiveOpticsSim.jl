include(joinpath(@__DIR__, "common.jl"))

import AdaptiveOpticsCalibration

const AOCMisregistration = AdaptiveOpticsCalibration.Misregistration

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

function estimate_raw_misregistration_offsets(meta::Calibration.MetaSensitivity,
    observed_response::AbstractMatrix)
    specification = AOCMisregistration.MetaSensitivityEstimateSpecification(
        meta.D0,
        meta.J,
        collect(meta.field_order),
        collect(meta.field_units),
    )
    plan = AdaptiveOpticsCalibration.prepare(
        AOCMisregistration.MetaSensitivityEstimate(), specification)
    result = AdaptiveOpticsCalibration.process(plan,
        AOCMisregistration.MetaSensitivityEstimateInputs(observed_response))
    return AOCMisregistration.parameter_offsets(result)
end

function apply_legacy_misregistration_offsets(meta::Calibration.MetaSensitivity,
    raw_offsets::AbstractVector; zero_point=Misregistration(T=eltype(raw_offsets)),
    gain::Real=1.0, precision::Int=4)
    estimate = zero_point
    for (index, field) in enumerate(meta.field_order)
        legacy_offset = round(gain * raw_offsets[index]; digits=precision)
        estimate = Calibration.update_misregistration(
            estimate,
            field,
            Calibration.misregistration_component(estimate, field) + legacy_offset,
        )
    end
    return estimate
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

    meta = Calibration.MetaSensitivity(reference, sensitivity_columns, epsilon, fields)
    injected = Misregistration(shift_x=5e-4, shift_y=-5e-4, T=Float64)
    dm_in = DeformableMirror(tel; n_act=dm.params.n_act,
        influence_model=influence_model(dm), misregistration=injected)
    response_in = geometric_wavefront_response_matrix(tel, dm_in, sensor, commands)
    raw_offsets = estimate_raw_misregistration_offsets(meta, response_in)
    estimate = apply_legacy_misregistration_offsets(meta, raw_offsets;
        zero_point=Misregistration(T=Float64), gain=1.0, precision=4)

    # A caller that needs another linearization explicitly reacquires `meta`
    # around `estimate`, measures its next response, then repeats the two calls
    # above. Neither AOS nor AOC retains zero-point iteration state.

    @info "Meta-sensitivity geometric-reference tutorial complete" shift_x=estimate.shift_x shift_y=estimate.shift_y
    return (
        injected=injected,
        estimate=estimate,
        raw_offsets=raw_offsets,
        response=reference,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

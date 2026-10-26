"""
    SecondOrderResponsePlan(; natural_frequency_hz, damping_ratio, sample_period_s,
                             T=Float32)

Run-immutable exact zero-order-hold discretization of
`y″ + 2ζωₙ y′ + ωₙ² y = ωₙ² u`, with `ωₙ = 2π * natural_frequency_hz`.
Frequency, damping ratio and sample period must be finite and positive. Input and
output have the same units; velocity has those units per second. Coefficients
are calculated in Float64, converted to `T` (Float32 or Float64), and rejected if
the converted homogeneous transition is not strictly stable. No delay,
saturation, instrument constants, or measured-response claim is implicit.
"""
struct SecondOrderResponsePlan{T<:Union{Float32,Float64}}
    natural_frequency_hz::Float64
    damping_ratio::Float64
    sample_period_s::Float64
    coefficients::NTuple{6,T}
end

_positive_response_parameter(value::Real) = isfinite(value) && value > 0
_positive_response_parameter(::Bool) = false

function SecondOrderResponsePlan(; natural_frequency_hz::Real,
        damping_ratio::Real, sample_period_s::Real, T::Type=Float32)
    T in (Float32, Float64) || throw(InvalidConfiguration(
        "second-order response precision must be Float32 or Float64"))
    parameters = (natural_frequency_hz, damping_ratio, sample_period_s)
    all(_positive_response_parameter, parameters) ||
        throw(InvalidConfiguration("second-order response parameters must be finite and positive"))
    frequency, damping, dt = Float64.(parameters)
    all(isfinite, (frequency, damping, dt)) || throw(InvalidConfiguration(
        "second-order response parameters must be representable in Float64"))
    omega = 2π * frequency
    # Augmented state [position, velocity, held input]; this also handles
    # critical damping without a branch or subtraction of nearly equal roots.
    generator = [0.0 dt 0.0; -dt*omega^2 -2dt*damping*omega dt*omega^2; 0.0 0.0 0.0]
    all(isfinite, generator) || throw(InvalidConfiguration(
        "second-order response generator exceeds Float64 range"))
    transition = exp(generator)
    coefficients = T.((transition[1,1], transition[1,2], transition[2,1],
        transition[2,2], transition[1,3], transition[2,3]))
    all(isfinite, coefficients) || throw(InvalidConfiguration(
        "second-order response coefficients exceed the execution precision"))
    a, b, c, d, _, _ = coefficients
    poles = eigvals([Float64(a) Float64(b); Float64(c) Float64(d)])
    all(p -> isfinite(p) && abs(p) < 1, poles) || throw(InvalidConfiguration(
        "second-order response transition must remain strictly stable at execution precision"))
    return SecondOrderResponsePlan{T}(frequency, damping, dt, coefficients)
end

"""Persistent position and velocity of independent response coordinates."""
struct SecondOrderResponseState{A}
    position::A
    velocity::A
end

_response_state_storage(::ScalarCPUStyle, input::AbstractVector{T}) where {T} =
    FixedSizeVectorDefault{T}(undef, length(input))
_response_state_storage(::AcceleratorStyle, input::AbstractVector) = similar(input)

"""
    PreparedSecondOrderResponse

Single-writer owner of exact input/output vector bindings, a response plan and
persistent scientific state. Construct with `prepare_second_order_response`.
"""
struct PreparedSecondOrderResponse{P,S,I,O}
    plan::P
    state::S
    input::I
    output::O
end

"""
    prepare_second_order_response(input, plan; output=similar(input))

Bind nonempty one-based Float32 or Float64 vectors on one compute device.
Input and output must have equal lengths, match the plan precision, and must
not alias. Position and velocity start at zero in independent persistent
storage. Mutate input values only between completed response advances.
"""
function prepare_second_order_response(input::AbstractVector{T},
        plan::SecondOrderResponsePlan{T}; output::AbstractVector{T}=similar(input)) where {T}
    Base.require_one_based_indexing(input, output)
    !isempty(input) && size(input) == size(output) || throw(DimensionMismatchError(
        "second-order response input/output must have equal nonzero lengths"))
    compute_device(input) == compute_device(output) || throw(InvalidConfiguration(
        "second-order response input/output must use the same compute device"))
    Base.mightalias(input, output) && throw(InvalidConfiguration(
        "second-order response input/output must not alias"))
    style = execution_style(input)
    position, velocity = _response_state_storage(style, input), _response_state_storage(style, input)
    owner = PreparedSecondOrderResponse(plan,
        SecondOrderResponseState(position, velocity), input, output)
    reset_response!(owner)
    return owner
end

@inline function _response_step(position, velocity, input, coefficients)
    a, b, c, d, g, h = coefficients
    return (a*position + b*velocity + g*input,
        c*position + d*velocity + h*input)
end

function _advance_response!(::ScalarCPUStyle, owner::PreparedSecondOrderResponse)
    position, velocity = owner.state.position, owner.state.velocity
    for i in eachindex(position, velocity, owner.input, owner.output)
        y, v = _response_step(position[i], velocity[i], owner.input[i], owner.plan.coefficients)
        position[i], velocity[i], owner.output[i] = y, v, y
    end
    return owner.output
end

@kernel function _second_order_response_kernel!(position, velocity, output, input, coefficients)
    i = @index(Global, Linear)
    y, v = _response_step(position[i], velocity[i], input[i], coefficients)
    position[i], velocity[i], output[i] = y, v, y
end

function _advance_response!(style::AcceleratorStyle, owner::PreparedSecondOrderResponse)
    launch_kernel!(style, _second_order_response_kernel!, owner.state.position,
        owner.state.velocity, owner.output, owner.input, owner.plan.coefficients;
        ndrange=length(owner.input))
    return owner.output
end

@inline _enqueue_response!(::ScalarCPUStyle, owner::PreparedSecondOrderResponse) =
    _advance_response!(ScalarCPUStyle(), owner)

function _enqueue_response!(style::AcceleratorStyle, owner::PreparedSecondOrderResponse)
    launch_kernel_async!(style, _second_order_response_kernel!, owner.state.position,
        owner.state.velocity, owner.output, owner.input, owner.plan.coefficients;
        ndrange=length(owner.input))
    return owner.output
end

"""
    advance_response!(prepared)

Hold the bound input constant for one plan sample period and return the bound
output at the end of that interval. Update persistent position and velocity
once; no input history or additional pure delay is inferred.
"""
advance_response!(owner::PreparedSecondOrderResponse) =
    _advance_response!(execution_style(owner.input), owner)

"""
    reset_response!(prepared)

Restore zero position, velocity and output without replacing their storage.
Serialize reset with response execution.
"""
function reset_response!(owner::PreparedSecondOrderResponse)
    fill!(owner.state.position, zero(eltype(owner.input)))
    fill!(owner.state.velocity, zero(eltype(owner.input)))
    fill!(owner.output, zero(eltype(owner.input)))
    return owner
end

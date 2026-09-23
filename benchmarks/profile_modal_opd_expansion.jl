"""Profile the retained plant-graph modal OPD executor on CPU, CUDA, or AMDGPU."""

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Calibration: ModalOPDExpansionPlan, combine_basis!
import AdaptiveOpticsCalibration
using Profile
using Statistics

if length(ARGS) == 1 && ARGS[1] == "cuda"
    import CUDA
elseif length(ARGS) == 1 && ARGS[1] == "amdgpu"
    import AMDGPU
end

function complete_call_ns!(output, plan, coefficients, synchronize!)
    started = time_ns()
    combine_basis!(output, plan, coefficients)
    synchronize!()
    return time_ns() - started
end

function steady_call_bytes(output, plan, coefficients)
    combine_basis!(output, plan, coefficients)
    return @allocated combine_basis!(output, plan, coefficients)
end

function run_case(basis, pupil_support, coefficients, synchronize!, target;
    warmup::Int=16, samples::Int=50)
    rows, columns, modes = size(basis)
    output = similar(basis, eltype(basis), rows, columns)
    plan = ModalOPDExpansionPlan(basis, pupil_support)
    for _ in 1:warmup
        combine_basis!(output, plan, coefficients)
    end
    synchronize!()

    elapsed_ns = Vector{Int}(undef, samples)
    for sample in eachindex(elapsed_ns)
        elapsed_ns[sample] = complete_call_ns!(output, plan, coefficients, synchronize!)
    end
    bytes = steady_call_bytes(output, plan, coefficients)
    synchronize!()

    support_host = Array(pupil_support)
    basis_host = Array(basis)
    coefficients_host = Array(coefficients)
    modal = AdaptiveOpticsCalibration.ModalBases
    reference_plan = AdaptiveOpticsCalibration.prepare(
        modal.ModalOPDExpansion(),
        modal.ModalOPDExpansionSpecification(
            rows, columns, modes, support_host, eltype(basis_host)),
    )
    expected = AdaptiveOpticsCalibration.process(
        reference_plan,
        modal.ModalOPDExpansionInputs(basis_host, coefficients_host),
    ).opd
    actual = Array(output)
    maximum(abs, actual .- expected) <= 1f-5 || error(
        "prepared plant-graph OPD expansion disagrees with AOC reference",
    )

    println("target=", target, " julia=", VERSION,
        " dimensions=", (rows, columns, modes), " warmup=", warmup,
        " samples=", samples, " host_bytes_per_call=", bytes,
        " p50_ns=", round(Int, median(elapsed_ns)),
        " p99_ns=", round(Int, quantile(elapsed_ns, 0.99)),
        " max_reference_error=", maximum(abs, actual .- expected))

    if target == "cpu"
        Profile.clear()
        Profile.@profile for _ in 1:100_000
            combine_basis!(output, plan, coefficients)
        end
        println("cpu_profile_calls=100000")
        Profile.print(format=:flat, sortedby=:count, maxdepth=12)
    elseif target == "cuda" && get(ENV, "AOS_MODAL_PROFILE_TRACE", "0") == "1"
        CUDA.Profile.profile_externally() do
            for _ in 1:20
                combine_basis!(output, plan, coefficients)
                synchronize!()
            end
        end
    end
    return nothing
end

function main(target::String)
    rows = 64
    columns = 64
    modes = 8
    basis = reshape(sin.(Float32.(1:(rows * columns * modes))) .* 0.01f0,
        rows, columns, modes)
    support = trues(rows, columns)
    support[1, :] .= false
    coefficients = Float32.(range(-0.2, 0.3; length=modes))

    if target == "cpu"
        return run_case(basis, support, coefficients, () -> nothing, target)
    elseif target == "cuda"
        CUDA.functional() || error("CUDA device is unavailable")
        return run_case(CUDA.CuArray(basis), CUDA.CuArray(support),
            CUDA.CuArray(coefficients), CUDA.synchronize, target)
    elseif target == "amdgpu"
        AMDGPU.functional() || error("AMDGPU device is unavailable")
        return run_case(AMDGPU.ROCArray(basis), AMDGPU.ROCArray(support),
            AMDGPU.ROCArray(coefficients), AMDGPU.synchronize, target)
    end
    error("expected target cpu, cuda, or amdgpu")
end

length(ARGS) == 1 || error("usage: julia profile_modal_opd_expansion.jl cpu|cuda|amdgpu")
main(ARGS[1])

using AdaptiveOpticsSim
using DelimitedFiles
using FFTW
using KernelAbstractions
using LinearAlgebra
using Random
using TOML

for name in names(AdaptiveOpticsSim; all=true)
    spelling = String(name)
    if Base.isidentifier(spelling) && !startswith(spelling, "#") &&
        !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

test_root = normpath(joinpath(@__DIR__, "..", "test"))
include(joinpath(test_root, "reference_harness.jl"))
include(joinpath(test_root, "gate0_characterization_harness.jl"))

root = isempty(ARGS) ? default_gate0_reference_root() : abspath(first(ARGS))
bundle = load_gate0_reference_bundle(root)
for case in bundle.cases
    case.baseline === :julia_gate0 || error(
        "refusing to generate non-Gate-0 case $(case.id)")
    actual = compute_gate0_actual(case)
    size(actual) == case.shape || error(
        "computed shape $(size(actual)) does not match $(case.shape) for $(case.id)")
    mkpath(dirname(case.data_path))
    write_reference_array(case.data_path, actual)
    println("wrote ", case.id, " -> ", case.data_path)
end

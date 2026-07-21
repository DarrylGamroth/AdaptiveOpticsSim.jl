import KernelAbstractions

# Shared test-only backend style for suites that compare scalar CPU and
# KernelAbstractions CPU implementations. Suites declare this fixture
# explicitly rather than relying on `ka_cpu_matrix.jl` to run first.
const KA_CPU_STYLE =
    AdaptiveOpticsSim.AcceleratorStyle(KernelAbstractions.CPU())

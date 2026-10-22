# Same physical model, old versus selected AOS batch policy.
# AOS_BATCH_VALIDATION_POLICY=baseline|selected (default selected).
const physical_batch_policy=get(ENV,"AOS_BATCH_VALIDATION_POLICY","selected")
physical_batch_policy in ("baseline","selected") || error("invalid batch validation policy")
if physical_batch_policy == "baseline"
    include("modulation_batch_override.jl")
    MODULATION_BATCH[]=0
end
include("validate_modulation_batches_copper_cuda.jl")
if abspath(PROGRAM_FILE) == @__FILE__
    batch_copper_main(ARGS)
end

include("pupil_modulation_storage_override.jl")
include("profile_physical_pyramid_cuda.jl")
if abspath(PROGRAM_FILE) == @__FILE__
    physical_profile_main(ARGS)
end

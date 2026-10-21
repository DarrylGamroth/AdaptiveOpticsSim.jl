# Profile one warmed pupil node in a complete captured plant step.
# AOS_SH_OCCUPANCY_GROUP=default or 44x8 selects the launch. Use ncu with
# --profile-from-start off --graph-profiling node and a pupil kernel filter.
import CUDA
include("benchmark_sh_pupil_workgroups_cuda.jl")

function profile_sh_pupil_occupancy()
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    selection = get(ENV, "AOS_SH_OCCUPANCY_GROUP", "default")
    selection in ("default", "44x8") || error("unknown occupancy workgroup")
    group = selection == "default" ? nothing : (44, 8, 1, 1)
    command, coordinates = _pdm_command(), _actuator_coordinates()
    opd = zeros(Float32, _PUPIL_RESOLUTION, _PUPIL_RESOLUTION)
    cpu = prepare_algorithm_graph(_plant_definition(command, opd, coordinates;
        name=:occupancy_cpu))
    _complete_step!(cpu)
    reference = copy(graph_output(cpu, Val(:frame)))
    gpu_command = CUDA.CuArray(command)
    definition = _plant_definition(gpu_command, CUDA.CuArray(opd),
        CUDA.CuArray(coordinates); name=:occupancy_cuda)
    fixture = sh_fixture(definition, Backends.compute_device(gpu_command),
        group, reference)
    println((; selection, geometry=sh_native_geometry(fixture.kernel.graph)))
    # All physical work is deterministic and the selected pupil kernel is
    # idempotent, so counter replay cannot consume a random or feedback history.
    CUDA.@profile external=true begin
        _complete_step!(fixture.graph)
    end
    _assert_reference_frame!(reference, fixture.graph, "post occupancy profile")
    println("occupancy_cpu_oracle=passed")
    return nothing
end

profile_sh_pupil_occupancy()

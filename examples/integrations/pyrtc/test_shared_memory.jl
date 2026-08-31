using Test

include(joinpath(@__DIR__, "pyrtc_shared_memory.jl"))
using .PyRTCSharedMemory

const PEER_PREFIX = "AOS_PYRTC_PEER "

function pyrtc_test_root()
    root = get(ENV, "PYRTC_ROOT", "")
    isdir(joinpath(root, "pyRTC")) || error(
        "set PYRTC_ROOT to a pyRTC checkout before running SHM integration tests",
    )
    return abspath(root)
end

function pyrtc_test_python()
    executable = get(
        ENV,
        "PYRTC_PYTHON",
        get(
            ENV,
            "JULIA_PYTHONCALL_EXE",
            something(Sys.which("python3"), ""),
        ),
    )
    isfile(executable) || error(
        "set PYRTC_PYTHON to the Python interpreter containing pyRTC dependencies",
    )
    return abspath(executable)
end

function peer_command(mode::Symbol, name::String, shape::Tuple)
    arguments = String[
        pyrtc_test_python(),
        joinpath(@__DIR__, "pyrtc_shm_peer.py"),
        String(mode),
        name,
    ]
    append!(arguments, string.(shape))
    append!(arguments, ("--pyrtc-root", pyrtc_test_root()))
    return Cmd(arguments)
end

function read_peer_control(process::Base.Process)
    while !eof(process)
        line = readline(process)
        startswith(line, PEER_PREFIX) || continue
        return chop(line; head=length(PEER_PREFIX), tail=0)
    end
    error("pyRTC interoperability peer exited without a control response")
end

function test_values(::Type{Float32}, shape::Tuple{Int})
    return Float32[index - 1 + 0.25 for index in 1:shape[1]]
end

function test_values(::Type{Float32}, shape::Tuple{Int,Int})
    return Float32[
        (row - 1) * 100 + (column - 1) + 0.25
        for row in 1:shape[1], column in 1:shape[2]
    ]
end

function close_and_unlink!(stream::PyRTCStream)
    close(stream)
    unlink!(stream)
    return nothing
end

function test_julia_producer_python_consumer(shape::Tuple)
    name = "aos_pyrtc_julia_producer_$(getpid())_$(length(shape))d"
    stream = create_stream(name, Float32, shape)
    process = nothing
    try
        process = open(peer_command(:consume, name, shape), "r+")
        @test read_peer_control(process) == "READY"
        @test publish!(stream, test_values(Float32, shape)) == UInt64(1)
        @test read_peer_control(process) == "DONE"
        wait(process)
        @test success(process)
    finally
        if !isnothing(process) && process_running(process)
            kill(process)
            wait(process)
        end
        close_and_unlink!(stream)
    end
    return nothing
end

function test_python_producer_julia_consumer(shape::Tuple)
    name = "aos_pyrtc_python_producer_$(getpid())_$(length(shape))d"
    process = open(peer_command(:produce, name, shape), "r+")
    stream = nothing
    try
        @test read_peer_control(process) == "READY"
        stream = open_stream(name, Float32, shape)
        write(process, "PUBLISH\n")
        flush(process)
        @test read_peer_control(process) == "PUBLISHED"
        actual = zeros(Float32, shape)
        @test read_next!(actual, stream; timeout=2.0) == UInt64(1)
        @test actual == test_values(Float32, shape)
        write(process, "STOP\n")
        flush(process)
        @test read_peer_control(process) == "DONE"
        wait(process)
        @test success(process)
    finally
        !isnothing(stream) && close(stream)
        if process_running(process)
            kill(process)
            wait(process)
        end
    end
    return nothing
end

@testset "native Julia and pyRTC ImageSHM interoperate" begin
    for shape in ((7,), (3, 5))
        test_julia_producer_python_consumer(shape)
        test_python_producer_julia_consumer(shape)
    end
end

@testset "native Julia pyRTC steady-state copies do not allocate" begin
    name = "aos_pyrtc_allocations_$(getpid())"
    owner = create_stream(name, Float32, (3, 5))
    peer = open_stream(name, Float32, (3, 5))
    source = test_values(Float32, (3, 5))
    destination = zeros(Float32, 3, 5)
    try
        publish!(owner, source)
        read_next!(destination, peer)

        publication_allocations = @allocated publish!(owner, source)
        read_allocations = @allocated read_next!(destination, peer)

        @test publication_allocations == 0
        @test read_allocations == 0
        @test destination == source
    finally
        close(peer)
        close_and_unlink!(owner)
    end
end

@testset "pyRTC SHM admission and ownership" begin
    name = "aos_pyrtc_admission_$(getpid())"
    owner = create_stream(name, Float32, (3, 5))
    peer = nothing
    try
        @test stream_shape(owner) == (3, 5)
        @test stream_count(owner) == UInt64(0)
        @test_throws PyRTCSharedMemoryError create_stream(
            name,
            Float32,
            (3, 5),
        )
        @test_throws PyRTCSharedMemoryError open_stream(
            name,
            Float64,
            (3, 5),
        )
        @test_throws PyRTCSharedMemoryError open_stream(
            name,
            Float32,
            (5, 3),
        )
        peer = open_stream(name, Float32, (3, 5))
        @test_throws PyRTCSharedMemoryError unlink!(peer)
        destination = zeros(Float32, 3, 5)
        @test_throws PyRTCSharedMemoryError read_next!(
            destination,
            peer;
            timeout=0.01,
            poll_interval=0.001,
        )
        @test_throws DimensionMismatch publish!(owner, zeros(Float32, 5, 3))
        @test_throws DimensionMismatch read_next!(
            zeros(Float32, 5, 3),
            peer,
        )
    finally
        !isnothing(peer) && close(peer)
        close_and_unlink!(owner)
    end

    @test_throws PyRTCSharedMemoryError create_stream(
        "invalid/name",
        Float32,
        (1,),
    )
    @test_throws PyRTCSharedMemoryError create_stream(
        "aos_pyrtc_zero_$(getpid())",
        Float32,
        (0,),
    )
    @test_throws PyRTCSharedMemoryError create_stream(
        "aos_pyrtc_rank_$(getpid())",
        Float32,
        (1, 1, 1),
    )

    collision_name = "aos_pyrtc_collision_$(getpid())"
    metadata_owner = create_stream(
        collision_name * "_meta",
        UInt8,
        (1,),
    )
    try
        @test_throws PyRTCSharedMemoryError create_stream(
            collision_name,
            Float32,
            (4,),
        )
        @test publish!(metadata_owner, UInt8[0x2a]) == UInt64(1)
        @test !ispath(joinpath("/dev/shm", collision_name))
    finally
        close_and_unlink!(metadata_owner)
    end
end

@testset "pyRTC SHM rejects malformed metadata and payload extent" begin
    metadata_name = "aos_pyrtc_bad_metadata_$(getpid())"
    metadata_owner = create_stream(metadata_name, Float32, (4,))
    try
        metadata_owner.metadata[4] = 99.0
        @test_throws PyRTCSharedMemoryError open_stream(
            metadata_name,
            Float32,
            (4,),
        )
    finally
        close_and_unlink!(metadata_owner)
    end

    extent_name = "aos_pyrtc_bad_extent_$(getpid())"
    extent_owner = create_stream(extent_name, Float32, (4,))
    close(extent_owner)
    open(joinpath("/dev/shm", extent_name), "r+") do io
        truncate(io, 3 * sizeof(Float32))
    end
    try
        @test_throws PyRTCSharedMemoryError open_stream(
            extent_name,
            Float32,
            (4,),
        )
    finally
        unlink!(extent_owner)
    end
end

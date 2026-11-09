# Optional selector, with CUDA and this checkout in the selected environment:
# julia --startup-file=no --project=<environment> test/cuda_spectral_synthesis.jl
# Add --benchmark for a captured-only comparison at fitted SPIDERS cardinalities.
using Test, AdaptiveOpticsSim, CUDA, LinearAlgebra, Random, Statistics
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends

CUDA.functional() || error("CUDA hardware is required for spectral synthesis tests")
CUDA.allowscalar(false)
BLAS.set_num_threads(1)

function cuda_spectral_node_types(graph)
    driver = CUDA.CUDACore
    count = Ref{Csize_t}(0)
    driver.cuGraphGetNodes(graph, C_NULL, count)
    nodes = Vector{driver.CUgraphNode}(undef, count[])
    driver.cuGraphGetNodes(graph, nodes, count)
    return map(nodes) do node
        kind = Ref{driver.CUgraphNodeType}()
        driver.cuGraphNodeGetType(node, kind)
        return kind[]
    end
end

function cuda_spectral_fixture(::Type{T}; summation=:coherent_modes,
        uniform=false, strided=false, shape=(7, 9), n_channel=5,
        n_anchor=3, n_path=4, cancellation=false, termwise=false) where {T}
    rng = MersenneTwister(9021)
    host_envelopes = [randn(rng, Complex{T}, shape) for _ in 1:(n_anchor * n_path)]
    basis = randn(rng, T, n_channel, n_anchor)
    coefficients = randn(rng, Complex{T}, n_channel, n_path)
    visibility = map(1:n_channel) do _
        factor = randn(rng, Complex{T}, n_path, n_path)
        factor ./= sqrt.(sum(abs2, factor; dims=2))
        factor * factor'
    end
    if termwise
        @assert n_path == 1 && n_anchor == 1 && uniform
        fill!(host_envelopes[1], Complex{T}(0.1, 0.2))
        fill!(basis, one(T))
        fill!(coefficients, one(Complex{T}))
        visibility = [ones(Complex{T}, 1, 1) for _ in 1:n_channel]
    elseif cancellation
        @assert n_path == 2 && n_anchor == 1
        fill!(host_envelopes[1], one(Complex{T}))
        fill!(host_envelopes[2], -one(Complex{T}))
        fill!(basis, one(T))
        fill!(coefficients, one(Complex{T}))
        # Pair arithmetic must retain this admitted negative PSD roundoff.
        visibility = [Complex{T}[1 1 + 2eps(T); 1 + 2eps(T) 1]
            for _ in 1:n_channel]
    end
    weights = uniform ? ones(T, n_channel) : T.(1:n_channel)
    wavelengths = collect(range(T(600e-9), T(800e-9); length=n_channel))
    channels = SpectralBundle(wavelengths, weights; T=T)
    response = uniform ? ones(T, n_channel) : collect(range(T(0.2), T(1.3);
        length=n_channel))
    !uniform && (response[1] = zero(T))
    cpu = Optics.prepare_spectral_synthesis(host_envelopes, basis, coefficients,
        visibility, channels; response, summation)
    staging = CuArray.(host_envelopes)
    device_envelopes = if strided
        map(host_envelopes) do envelope
            storage = CUDA.fill(Complex{T}(NaN, NaN), 2shape[1], 2shape[2])
            view(storage, 2:2:2shape[1], 2:2:2shape[2])
        end
    else
        CuArray.(host_envelopes)
    end
    output_parent = CUDA.fill(T(NaN), strided ? (2shape[1], 2shape[2]) : shape)
    output = strided ? view(output_parent, 2:2:2shape[1], 2:2:2shape[2]) : output_parent
    for (destination, source) in zip(device_envelopes, staging)
        destination .= source
    end
    gpu = Optics.prepare_spectral_synthesis(device_envelopes, basis, coefficients,
        visibility, channels; response, summation, output)
    # Both executors must use their frozen snapshots, including resident copies.
    fill!(basis, zero(T)); fill!(coefficients, zero(Complex{T}))
    fill!(response, zero(T)); foreach(matrix -> fill!(matrix, zero(Complex{T})), visibility)
    return (; cpu, gpu, host_envelopes, staging, output_parent, strided)
end

function cuda_spectral_copy_inputs!(fixture)
    for (host, staging, bound) in zip(fixture.host_envelopes, fixture.staging,
            fixture.gpu.envelopes)
        copyto!(staging, host)
        bound .= staging
    end
    return nothing
end

function cuda_spectral_poison_scratch!(owner)
    fill!(owner.products.output, NaN)
    nan = eltype(owner.workspace.coherent_field)(NaN, NaN)
    foreach(field -> fill!(field, nan), owner.workspace.channel_fields)
    fill!(owner.workspace.coherent_field, nan)
    fill!(owner.workspace.packed_envelopes, nan)
    fill!(owner.workspace.packed_channel_fields, nan)
    return nothing
end

function cuda_spectral_capture_pair(owner, context)
    # Streamed work is confined to warmup and numerical gates; measured work
    # below always replays captured device graphs on this retained stream.
    Optics._synthesize_spectrum_reference!(owner)
    Optics.synthesize_spectrum!(owner)
    CUDA.synchronize(context.stream)
    reference = Backends._capture_prepared_device_graph(context) do
        Optics._synthesize_spectrum_reference!(owner)
        nothing
    end
    fused = Backends._capture_prepared_device_graph(context) do
        Optics.synthesize_spectrum!(owner)
        nothing
    end
    kinds = cuda_spectral_node_types(fused.graph)
    @test length(kinds) == length(owner.envelopes) + 2
    @test all(==(CUDA.CUDACore.CU_GRAPH_NODE_TYPE_KERNEL), kinds)
    reference_nodes = length(cuda_spectral_node_types(reference.graph))
    @test reference_nodes > length(kinds)
    println((; precision=eltype(owner.products.output), summation=owner.plan.summation,
        channels=length(owner.plan.channels), anchors=owner.plan.n_anchor,
        paths=owner.plan.n_path, reference_nodes, fused_nodes=length(kinds)))
    flush(stdout)
    return (; reference, fused)
end

cuda_spectral_host_allocations(executable) = @allocated CUDA.launch(executable)
cuda_spectral_device_allocations(executable) = CUDA.@allocated CUDA.launch(executable)

function compare_cuda_spectral_capture(::Type{T}; cancellation=false, termwise=false,
        kwargs...) where {T}
    fixture = cuda_spectral_fixture(T; cancellation, termwise, kwargs...)
    owner = fixture.gpu
    # The context seam accepts owning CUDA storage; the bound output may be a view.
    context = Backends._prepare_device_execution_context(fixture.output_parent)
    retained = (owner.plan.resident_calibration, owner.workspace.packed_envelopes,
        owner.workspace.packed_channel_fields, owner.products.output, owner.envelopes)
    original_inputs = copy.(fixture.host_envelopes)
    Backends._with_prepared_device_execution_context(context) do
        cuda_spectral_copy_inputs!(fixture)
        graphs = cuda_spectral_capture_pair(owner, context)
        first_image = nothing
        for frame in 1:3
            if frame == 2
                changed_input = fixture.host_envelopes[min(2, length(fixture.host_envelopes))]
                changed_input .*= Complex{T}(0.3, -0.8)
            elseif frame == 3
                # Synthesis has no persistent scientific state or reset API.
                # Restore inputs and poison every scratch/output to check replay
                # overwrites all prior-frame values without recapturing.
                foreach(pair -> copyto!(pair...), zip(fixture.host_envelopes, original_inputs))
                cuda_spectral_poison_scratch!(owner)
            end
            cuda_spectral_copy_inputs!(fixture)
            expected = copy(Optics.synthesize_spectrum!(fixture.cpu))
            Optics._synthesize_spectrum_reference!(owner)
            CUDA.synchronize(context.stream)
            original_gpu = Array(owner.products.output)
            Optics.synthesize_spectrum!(owner)
            CUDA.synchronize(context.stream)
            streamed_gpu = Array(owner.products.output)
            images = map((graphs.reference, graphs.fused)) do graph
                fill!(owner.products.output, T(NaN))
                Backends._launch_prepared_device_graph!(graph, context)
                CUDA.synchronize(context.stream)
                Array(owner.products.output)
            end
            for actual in (original_gpu, streamed_gpu, images...)
                @test all(isfinite, actual)
                @test actual ≈ expected rtol=128eps(T) atol=128eps(T)
                @test actual ≈ original_gpu rtol=128eps(T) atol=128eps(T)
                if termwise
                    @test actual == original_gpu
                end
            end
            if frame == 1
                first_image = copy(last(images))
                if cancellation
                    if owner.plan.summation === :path_pairs
                        @test all(<(zero(T)), last(images))
                    else
                        @test all(>=(zero(T)), last(images))
                        @test maximum(last(images)) <= 64(eps(Float64) + Float64(eps(T))^2)
                    end
                end
            elseif frame == 2
                @test !isapprox(last(images), first_image; rtol=1e-4, atol=128eps(T))
            else
                @test last(images) ≈ first_image rtol=128eps(T) atol=128eps(T)
            end
            @test (owner.plan.resident_calibration, owner.workspace.packed_envelopes,
                owner.workspace.packed_channel_fields, owner.products.output,
                owner.envelopes) === retained
        end
        for graph in (graphs.reference, graphs.fused)
            CUDA.launch(graph.executable)
            cuda_spectral_host_allocations(graph.executable)
            cuda_spectral_device_allocations(graph.executable)
            @test cuda_spectral_host_allocations(graph.executable) == 0
            @test cuda_spectral_device_allocations(graph.executable) == 0
            CUDA.synchronize(context.stream)
        end
        if fixture.strided
            @test all(isnan, Array(fixture.output_parent)[1:2:end, :])
            @test all(isnan, Array(fixture.output_parent)[:, 1:2:end])
        end
    end
    return nothing
end

@testset "CUDA spectral synthesis same-owner streamed and captured oracles" begin
    for T in (Float32, Float64), summation in (:coherent_modes, :path_pairs),
            uniform in (false, true), strided in (false, true)
        compare_cuda_spectral_capture(T; summation, uniform, strided)
    end
    for T in (Float32, Float64), summation in (:coherent_modes, :path_pairs)
        compare_cuda_spectral_capture(T; summation, uniform=true, n_channel=3,
            n_anchor=1, n_path=2, cancellation=true)
        compare_cuda_spectral_capture(T; summation, uniform=true, n_channel=3,
            n_anchor=1, n_path=1, termwise=true)
    end
end

function cuda_spectral_replay_seconds(graph, context; repetitions=31)
    for _ in 1:10
        Backends._launch_prepared_device_graph!(graph, context)
    end
    CUDA.synchronize(context.stream)
    return [@elapsed begin
        Backends._launch_prepared_device_graph!(graph, context)
        CUDA.synchronize(context.stream)
    end for _ in 1:repetitions]
end

if "--benchmark" in ARGS
    @testset "CUDA captured spectral synthesis comparison" begin
        for T in (Float32, Float64), summation in (:coherent_modes, :path_pairs)
            fixture = cuda_spectral_fixture(T; summation, uniform=true,
                shape=(127, 129), n_channel=8, n_anchor=3, n_path=13)
            context = Backends._prepare_device_execution_context(fixture.gpu.products.output)
            Backends._with_prepared_device_execution_context(context) do
                graphs = cuda_spectral_capture_pair(fixture.gpu, context)
                reference = cuda_spectral_replay_seconds(graphs.reference, context)
                fused = cuda_spectral_replay_seconds(graphs.fused, context)
                println((; T, summation, shape=(127, 129), channels=8, anchors=3,
                    paths=13, reference_nodes=length(cuda_spectral_node_types(graphs.reference.graph)),
                    fused_nodes=length(cuda_spectral_node_types(graphs.fused.graph)),
                    reference_median_ms=1e3median(reference), fused_median_ms=1e3median(fused),
                    reference_p95_ms=1e3quantile(reference, 0.95),
                    fused_p95_ms=1e3quantile(fused, 0.95)))
                flush(stdout)
            end
        end
    end
end

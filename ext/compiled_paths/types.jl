# The compiler's centred (y,x) arrays are the numerical contract. Only PROPER
# propagation crosses to its corner-origin storage, through prepared scratch.
@enum PathOperationKind begin
    PathMultiply
    PathNormalize
    PathTransferFunction
    PathFourierTransform
    PathSpectralMultiplier
    PathResample
    PathLiveScreen
end

struct PathOperation{T,A}
    kind::PathOperationKind
    array::A
    auxiliary::A
    post::A
    dx::T
    dz::T
    wavelength_m::T
    gouy::Complex{T}
    from_waist::Bool
    forward::Bool
    source::UInt8 # 1 = DM OPD; 2 = entrance OPD and log-amplitude attenuation
end

# Closed execution vocabulary. Geometry and cache indices remain values; only
# the finite transform-group arities specialize tuple arithmetic. This avoids
# one Julia type/method family per complete YAML prescription.
@enum PathExecutionKind begin
    ExecuteMultiply
    ExecuteNormalize
    ExecutePropagation
    ExecuteSpectralMultiplier
    ExecuteResample
    ExecuteLiveDM
    ExecuteLiveEntrance
    ExecutePropagationPair
    ExecuteGroup2
    ExecuteGroup4
    ExecuteGroup6
    ExecuteGroup8
    ExecuteGroup10
    ExecutePropagationResample
    ExecuteCancelledResample
    LoadCheckpoint
    StoreCheckpoint
end

struct PathExecutionStep
    kind::PathExecutionKind
    index::Int # operation position within its sequence
    mask::Int  # initial multiplier position for a transform group, or zero
    cache::Int # checkpoint slot; zero means the current field for a group
end

struct PathSequenceExecution
    first::Int
    last::Int
end

# Immutable execution decisions only. Exact arrays, FFT plans and live inputs
# are resolved through the prepared owner, never retained in this contract.
struct PathExecutionPlan{S,R}
    steps::S
    shared::R
    unshared::R
end

struct OpticalPathPlan{A,O}
    path_index::Int
    anchor_index::Int
    entrance::A
    prefix::O
    closed::O
    only::O
end

struct PathTransformPlan{A,V}
    source_index::Int
    shift_transfer::A
    phase_multiplier::V
    defocus_transfer::V
    ghost::V
end

struct PathPrefixCheckpoint
    index::Int
    cache_index::Int
end

struct PathPrefixPlan{C}
    cache_index::Int
    count::Int
    checkpoints::C
end

struct PathSharingPlan{S}
    prefix::S
    closed::S
    only::S
end

struct CompiledPathsPlan{T,D,P,F,A,S}
    target::D
    grid::Int
    image_size::Int
    crop_start::Int
    n_env::Int
    n_anchor::Int
    wavelengths::FixedSizeVectorDefault{T}
    paths::P
    shared_prefixes::S
    fast_paths::F
    resample_kernel::A
    main_phase::A
    pinhole_gain::FixedSizeVectorDefault{Complex{T}}
    pixel_mtf::A
    phantom_pad::Int
    pixel_pad::Int
    apply_pixel_mtf::Bool
end

struct PathOperationWorkspace{A,R,E}
    input::A
    temporary::A
    mapped::A
    opd::R
    electric_field::E
    identity_roll::Bool
end

struct PathFFTWorkspace{A,F,I}
    field::A
    scratch::A
    forward::F
    inverse::I
end

# A plane is one fast path at one spectral anchor. Transforms cover only y,x;
# the third axis is an independent batch, never an optical coordinate.
struct PathFilterBatch{A,F,I}
    field::A
    transfer::A
    forward::F
    inverse::I
end

struct PathFastWorkspace{B,D,S}
    image::B
    defocus::D
    defocus_slots::S
end

struct PathNormalizationWorkspace{P,S}
    partials::P
    power::S
end

struct PathResampleBatch{N,A,F,O,C}
    batch_size::Val{N}
    field::A
    temporary::A
    kernel::A
    forward::F
    outputs::O
    coefficients::C
end

struct PathBranchTransform{Terminal}
    index::Int
    mask::Int
end

struct PathBranchSegment
    last::Int
    closed_cache::Int
    only_cache::Int
    transforms::FixedSizeVectorDefault{Union{PathBranchTransform{false},PathBranchTransform{true}}}
end

struct PathBranchPair{S}
    first::Int
    segments::S
end

struct PathBranchBatch{A,V,F,B,P}
    field::A
    planes::V
    forward::F
    backward::B
    pairs::P
end

struct CompiledPathsWorkspace{A,W,C,O,F,V,N,M,L,B,P,E,K,R,Q}
    field::A
    prefix::A
    shared_prefixes::V
    wavefront::W
    context::C
    operations::O
    physical_fft::F
    image_fft::F
    phantom_fft::F
    crop::A
    resample_temporary::A
    closed_envelopes::V
    only_envelopes::V
    open_envelopes::V
    closed_sources::V
    only_sources::V
    normalization::N
    pixel_mtf::M
    mul_coefficients::L
    fast::B
    propagation::P
    edge_multiplier::E
    prepared_mul_coefficients::K
    image_resample::R
    branch_batch::Q
end

"""Canonical `(y,x)` compiled multi-path relative-intensity products, dimensionless intensity."""
struct CompiledPathsProducts{R}
    fringed::R
    unfringed::R
    reference::R
end

"""Prepared compiled optical path execution with exact input buffer bindings and products."""
struct PreparedCompiledPaths{P,W,I,O,C,S,E}
    plan::P
    workspace::W
    inputs::I
    products::O
    execution_context::C
    synthesis::S
    execution::E
end

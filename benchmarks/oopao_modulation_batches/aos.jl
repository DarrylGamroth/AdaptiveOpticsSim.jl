# Usage: julia --project=ENV aos.jl [REFERENCE_DIRECTORY]
# AOS_ORACLE_BACKEND=cpu|cuda; AOS_ORACLE_BATCH optional benchmark override.
using AdaptiveOpticsSim, DelimitedFiles, LinearAlgebra, TOML
using AdaptiveOpticsSim.Optics, AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Backends
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const W = AdaptiveOpticsSim.WavefrontSensors
const oracle_root = isempty(ARGS) ? get(ENV, "AOS_ORACLE_OUTPUT", @__DIR__) : only(ARGS)
const oracle_backend = get(ENV, "AOS_ORACLE_BACKEND", "cpu")
oracle_backend in ("cpu", "cuda") || error("AOS_ORACLE_BACKEND must be cpu or cuda")
const oracle_requested_batch = parse(Int, get(ENV, "AOS_ORACLE_BATCH", "0"))
if haskey(ENV, "AOS_ORACLE_BATCH")
    default_override = joinpath(dirname(dirname(pathof(AdaptiveOpticsSim))), "benchmarks", "modulation_batch_override.jl")
    include(get(ENV, "AOS_ORACLE_BATCH_OVERRIDE", default_override))
    MODULATION_BATCH[] = oracle_requested_batch
end
if oracle_backend == "cuda"
    @eval import CUDA
    CUDA.allowscalar(false)
end
const oracle_selector = oracle_backend == "cuda" ? CUDABackend() : CPUBackend()
const oracle_precision = get(ENV, "AOS_ORACLE_PRECISION", "Float64")
oracle_precision in ("Float32", "Float64") || error("precision must be Float32 or Float64")
const oracle_T = oracle_precision == "Float32" ? Float32 : Float64
const oracle_strategy = get(ENV, "AOS_ORACLE_STRATEGY", "pupil_tilt")
oracle_strategy in ("pupil_tilt", "shifted_mask") || error("invalid oracle strategy")
const precision_suffix = oracle_T === Float32 ? "_f32" : ""
const oracle_default_prefix = (oracle_backend == "cpu" && !haskey(ENV,"AOS_ORACLE_BATCH") ? "aos" : "aos_$(oracle_backend)_b$(oracle_requested_batch)") * precision_suffix
const oracle_prefix = get(ENV, "AOS_ORACLE_PREFIX", oracle_default_prefix)
for (n,p,m,points,sep,edge,b) in ((24,4,1.,8,4,0,2),(32,8,3.,20,8,0,1),(120,30,3.,32,2,1,1))
    name = "n$(n)_p$(p)_m$(Int(m))"
    selected = get(ENV, "AOS_ORACLE_CASES", "n24_p4_m1,n32_p8_m3,n120_p30_m3")
    name in split(selected, ',') || continue
    path = joinpath(oracle_root, name)
    tel = Telescope(resolution=n, diameter=8., central_obstruction=0., T=oracle_T, backend=oracle_selector)
    src = Source(band=:I, magnitude=0., T=oracle_T)
    pupil = PupilFunction(tel; T=oracle_T, backend=oracle_selector)
    wfs = PyramidWFS(tel; pupil_samples=p, modulation=m, modulation_points=points,
        n_pix_separation=sep, n_pix_edge=edge, binning=b,
        modulation_propagation_strategy=oracle_strategy == "pupil_tilt" ?
            W.PyramidPupilTiltStrategy() : W.PyramidShiftedMaskStrategy(), T=oracle_T)
    front = PyramidOpticalFrontEnd(wfs, src)
    rate = pyramid_rate_map(front, pupil)
    prepared = prepare_wfs_optics(front, pupil, rate)
    batch = W.pyramid_propagation_workspace(wfs).modulation_batch
    metadata = Dict("precision"=>string(oracle_T),"backend"=>oracle_backend,
        "actual_batch"=>hasproperty(batch,:batch_size) ? batch.batch_size : 1,
        "padded_side"=>size(W.pyramid_propagation_workspace(wfs).field,1),
        "points"=>points,"strategy"=>oracle_strategy)
    open(io->TOML.print(io,metadata;sorted=true),joinpath(path,oracle_prefix*"_parameters.toml"),"w")
    function frame(opd)
        copyto!(pupil.opd, oracle_T.(opd))
        form_wfs_optical_products!(rate, pupil, prepared)
        return Array(rate.values) .* .001
    end
    flat = frame(zeros(n,n))
    ramp = frame([5e-9*(i-1)-2e-9*(j-1) for i in 1:n,j in 1:n])
    basis = repeat(reshape(collect(range(-1.; step=2/n, length=n)),n,1),1,n)
    plus = frame(1e-9 .* basis)
    minus = frame(-1e-9 .* basis)
    derivative = (plus-minus)/2e-9
    for (label,a) in (("flat_frame",flat),("ramp_frame",ramp),("tip_frame_derivative",derivative))
        writedlm(joinpath(path,oracle_prefix*"_"*label*".txt"),a)
        ref = readdlm(joinpath(path,label*".txt"))
        flat_ref = readdlm(joinpath(path,"flat_frame.txt"))
        println(name," ",oracle_prefix," ",label," flux=",sum(a)," ref=",sum(ref),
            " relnorm=",norm(a-ref)/norm(ref),
            " normalized_relnorm=",norm(a/sum(flat)-ref/sum(flat_ref))/(norm(ref)/sum(flat_ref)))
    end
end

"""Archive the AOS implementation experiment and reduce native CUDA activity."""
from pathlib import Path
import collections,gzip,hashlib,json,shutil,statistics,subprocess,tomllib
base=Path(__file__).resolve().parent
repo=Path('/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-copper-review')
out=repo/'benchmarks/results/platform/2026-10-21-copper-implementation'
raw=out/'raw';raw.mkdir(parents=True,exist_ok=True)
for source in (base/'aos-profiles').iterdir():
    if source.suffix=='.log':
        (raw/(source.name+'.gz')).write_bytes(gzip.compress(source.read_bytes(),mtime=0))
    else:shutil.copy2(source,raw/source.name)
for name in ('native-trace.jl','run-profiles.py','profile-baseline-harness.jl','gpu-start.txt','gpu-end.txt','system.txt'):
    shutil.copy2(base/name,raw/name)
for name in ('Project.toml','Manifest.toml'):
    (raw/(name+'.gz')).write_bytes(gzip.compress((base/'profile-env'/name).read_bytes(),mtime=0))
shutil.copy2(base/'aos-source-review.md',out/'source-review.md')
records=[tomllib.loads((raw/f'run{i}.toml').read_text()) for i in (1,2,3)]
sha=hashlib.sha256((raw/'profile-baseline-harness.jl').read_bytes()).hexdigest()
for r in records:
    assert r['harness_sha256']==sha and r['repeats']==128
    assert r['captured']['host_bytes']==[0]*16
    for mode in ('stream','captured'):
        assert len(r[mode]['samples_ms'])==128
        assert statistics.median(r[mode]['samples_ms'])==r[mode]['median_ms']
summary={mode:dict(median_ms=statistics.median(r[mode]['median_ms'] for r in records),
    min_run_median_ms=min(r[mode]['median_ms'] for r in records),
    max_run_median_ms=max(r[mode]['median_ms'] for r in records),
    host_bytes=sorted({v for r in records for v in r[mode]['host_bytes']})) for mode in ('stream','captured')}
node_summary={mode:{name:statistics.median(r[mode][name]['device_event']['median_ms']*1000 for r in records)
    for name in records[0][mode]} for mode in ('nodes_stream','nodes_captured')}
activity=[];mode=None;device=False;busy={}
for line in (base/'aos-profiles/native-trace.log').read_text().splitlines():
    if line.startswith('Native complete-frame timeline:'):
        mode=line.rsplit(' ',1)[1];device=False
    if line.startswith('Device-side activity:'):
        device=True;busy[mode]=line
    if not device or not line.startswith('│'):continue
    columns=[v.strip() for v in line.split('│')[1:-1]]
    if len(columns)!=10 or not columns[0].isdigit():continue
    duration,unit=columns[2].split()
    activity.append(dict(mode=mode,correlation_id=int(columns[0]),start=columns[1],
        duration_us=float(duration)*{'ns':.001,'µs':1,'ms':1000,'s':1e6}[unit],
        threads=columns[3],blocks=columns[4],registers=columns[5],shared_memory=columns[6],
        size=columns[7],throughput=columns[8],name=columns[9]))
assert activity and set(busy)=={'stream','captured'}
sums=collections.defaultdict(float);counts=collections.Counter()
for row in activity:
    if row['mode']=='captured':
        name=row['name'].split('(')[0]
        sums[name]+=row['duration_us'];counts[name]+=1
kernels=[dict(name=name,calls=counts[name],mean_us_per_frame=value/4)
         for name,value in sorted(sums.items(),key=lambda item:-item[1])]
(out/'summary.json').write_text(json.dumps(dict(complete=summary,node_device_intervals_us=node_summary,
    native_trace_busy=busy,native_captured_kernels=kernels),indent=2)+'\n')
(out/'activity.json').write_text(json.dumps(activity,indent=2)+'\n')
lines=['# AOS implementation costs for the Copper CUDA workload','',
    '2026-10-21, RTX 3080, driver 615.71.10. This holds the production 480²/32-point '
    'prescription fixed. No AOS numerical algorithm, instrument setting or package API changed.','',
    '## Completed frame timings','',
    'Three fresh Julia processes, 32 warmups and 128 samples per mode. Medians below are '
    'medians of process medians; the range is between process medians. The boundary is '
    '`step_graph!` through prepared completion, including native five-layer atmosphere and noisy '
    'detector, with zero DM commands. It excludes HIL staging, transport and the controller.','',
    '| Execution | Median ms | Run-median range ms | Warm host bytes/frame |','|---|---:|---:|---:|']
for mode in ('stream','captured'):
    s=summary[mode];lines.append(f'| {mode} | {s["median_ms"]:.3f} | {s["min_run_median_ms"]:.3f}–{s["max_run_median_ms"]:.3f} | {s["host_bytes"]} |')
lines+=['','## Native kernel trace','',busy['captured']+'. This is a four-frame profiled window, '
    'not whole-application utilization, achieved occupancy, peak efficiency or a fixed-rate test. '
    'Profiler durations are diagnostic and are not the unprofiled latency baseline.','',
    '| Captured activity | Mean µs/frame |','|---|---:|']
for k in kernels[:8]:lines.append(f'| `{k["name"]}` | {k["mean_us_per_frame"]:.2f} |')
lines+=['','Both cuFFT rows include forward and backward transforms. Most of their time is in the '
    '32-plane inverse batch. The shifted-mask fill and intensity reduction are the next large costs. '
    'The intensity kernel reads a 256 MiB complex stack plus a 4 MiB prior intensity and writes '
    '4 MiB; its observed duration implies about 700 GB/s of algorithmic array traffic. This is a '
    'derived byte rate, not a hardware DRAM-counter measurement. Removing a small final image pass '
    'does not remove the dominant complex-stack traffic.','',
    '## Serialized node attribution','',
    '| Node | Stream device interval µs | Captured-protocol device interval µs |',
    '|---|---:|---:|']
for name in ('atmosphere','pdm','pupil_opd_composition','pwfs','detector'):
    lines.append(f'| {name} | {node_summary["nodes_stream"][name]:.2f} | {node_summary["nodes_captured"][name]:.2f} |')
lines+=['','These use separate disposable graphs and an event around direct node submission. '
    'They include device idle gaps between submissions and are measured with per-node synchronization; captured '
    'attribution calls the capture protocol but does not launch the native graph. They are not '
    'additive stage costs or native replay scheduling measurements. The archived host-complete '
    'node times additionally use CUDA event waiting, whose task wait policy differs from the '
    'prepared graph completion policy; do not substitute them for full-frame latency.','',
    '## Implementation findings','',
    '1. PWFS propagation is the main warmed target. Benchmark the fixed batch heuristic at unchanged '
    'pupil/modulation sampling before selecting another value. Preserve weights, FFT normalization, '
    'all points, captured reset/replay behavior and numerical parity. No optimum is established here.',
    '2. Preparation allocates dense pupil modulation phases even when the shifted-mask rate path '
    'uses compact axis factors. At 1440/256 this is 3.955 GiB device storage and an equally large '
    'temporary host array. At the production 480/32 setting it is 56.25 MiB each. The host temporary '
    'is copied once during preparation; it is not a per-frame transfer. Required phases can be '
    'generated on device and remain resident. Strategy-specific ownership can avoid the cube '
    'entirely for the rate path while preserving diagnostic/reference products. No memory-leak '
    'claim or implemented speedup is made.',
    '3. Binning, atmospheric interpolation and DM kernels are secondary measured costs. They may '
    'benefit from output reductions or fusion, but cannot explain most of this frame time. The '
    'executor already replays one native graph; an executor rewrite is not the first target.','',
    'See the [source review](source-review.md) for affected code, existing efficient behavior and '
    'validation obligations. Kernel identifiers, grids, registers and shared-memory attributes '
    'remain in [activity.json](activity.json) and the full trace. No occupancy counters or '
    'roofline experiment were collected.','',
    '## Reproduction and provenance','',
    'Use an isolated environment with local AOS/Copper sources. The exact benchmark commands and '
    'environment snapshots are archived under raw/. Julia 1.13.2, CUDA 6.4.2, KA 0.9.45, CPU 2 '
    'affinity, one compute/GC/BLAS/FFT thread, no interactive pool, shared host and uncontrolled '
    'GPU clocks. These are sequential service times, not an open-load latency contract.','',
    'Timed runs used the archived `profile-baseline-harness.jl` from AOS commit 041cb2a2. '
    'A subsequent harness-only change displays the CUDA 6.4 returned profiling object explicitly; '
    'it does not change timed paths or numerical source. The native trace uses its separately '
    'archived driver. All raw timestamps in activity.json are reductions of the human-readable '
    'trace, with its display rounding. No registered application environment was modified.','']
(out/'README.md').write_text('\n'.join(lines))
shutil.copy2(Path(__file__),raw/'archive-aos-profile.py')
provenance=dict(baseline_harness_sha256=sha,current_harness_sha256=hashlib.sha256((repo/'benchmarks/profile_copper_cuda.jl').read_bytes()).hexdigest(),
    numerical_source_tree=subprocess.check_output(['git','rev-parse','HEAD:src'],cwd=repo,text=True).strip(),
    extension_tree=subprocess.check_output(['git','rev-parse','HEAD:ext'],cwd=repo,text=True).strip(),
    raw_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in raw.iterdir() if p.is_file()})
(out/'provenance.json').write_text(json.dumps(provenance,indent=2,sort_keys=True)+'\n')
print('Archived AOS profile:',len(activity),'device activities and',len(provenance['raw_sha256']),'raw artifacts')

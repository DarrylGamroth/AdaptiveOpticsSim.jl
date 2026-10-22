"""Record shared-host conditions; never change clocks, power or other jobs."""
import argparse
import json
import signal
import subprocess
import time
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('output',type=Path)
a=p.parse_args();running=True

def stop(*_):
    global running
    running=False

signal.signal(signal.SIGTERM,stop)
signal.signal(signal.SIGINT,stop)
with a.output.open('x') as f:
    while running:
        gpu=subprocess.run(['nvidia-smi','--query-gpu=timestamp,pstate,utilization.gpu,utilization.memory,memory.used,temperature.gpu,power.draw,clocks.sm,clocks.mem','--format=csv,noheader,nounits'],text=True,capture_output=True)
        processes=subprocess.run(['ps','-eo','pid,psr,pcpu,comm','--sort=-pcpu'],text=True,capture_output=True)
        record=dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            gpu_csv=gpu.stdout.strip(),gpu_status=gpu.returncode,
            loadavg=Path('/proc/loadavg').read_text().strip(),
            procstat=Path('/proc/stat').read_text().splitlines()[:17],
            processes=processes.stdout.splitlines()[:20])
        f.write(json.dumps(record)+'\n');f.flush()
        time.sleep(2)

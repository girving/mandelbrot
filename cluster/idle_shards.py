#!/usr/bin/env python3
"""Run sharded escape_tree runs on idle GPUs only, one single-GPU job per shard, then merge each run.

  cluster/idle_shards.py RUNS.json

RUNS.json lists runs: [{"name": ..., "shards": N, "command": "./build/release/escape_tree --cuda ... ks",
"build64": false, "timeout": seconds per shard, "free_gpus": G, "max_active": M}, ...].  The optional
free_gpus and max_active override FREE_GPUS (GPUs that must stay unused for a shard to start) and cap that
run's concurrently active shards.  Shard results persist in /data/results/<name> on the
PVC (cluster/bench.sh's SHARDED mode), so rerunning this after an interruption only redoes missing shards.

The GPU queue has no priority classes and no preemption, so politeness is enforced here: a shard is submitted
only while the queue has no pending workloads, at least FREE_GPUS GPUs stay unused (SHORT_FREE_GPUS for runs
marked "short": true, which take a few minutes), and at most MAX_JOBS of our shard jobs are active.  Long jobs
cannot be preempted, so they start only when the cluster is far from full.  Once a run's shards are all done, a merge job (one GPU, a few minutes: it builds,
then only merges) prints the whole run's report; its log is the result.
"""
import json
import subprocess
import sys
import time

NAMESPACE, QUEUE, FREE_GPUS, SHORT_FREE_GPUS, MAX_JOBS, TOTAL_GPUS = 'research', 'gpu', 64, 16, 8, 248
REPO = 'girving/mandelbrot'


def kubectl(*args, input=None, check=True):
    r = subprocess.run(['kubectl', '-n', NAMESPACE, *args], input=input, capture_output=True, text=True, timeout=60)
    if check and r.returncode:
        raise RuntimeError(f'kubectl {" ".join(args)}: {r.stderr.strip()}')
    return r.stdout


def queue_state():
    """(pending workloads, GPUs in use) for the GPU cluster queue"""
    q = json.loads(subprocess.run(['kubectl', 'get', 'clusterqueue', QUEUE, '-o', 'json'], capture_output=True,
                                  text=True, timeout=60, check=True).stdout)
    used = 0
    for f in q['status'].get('flavorsUsage', []):
        for r in f['resources']:
            if r['name'] == 'nvidia.com/gpu':
                used += int(r['total'])
    return q['status'].get('pendingWorkloads', 0), used


def job_status(name):
    """'active', 'complete', 'failed', or None if absent"""
    r = subprocess.run(['kubectl', '-n', NAMESPACE, 'get', 'job', name, '-o', 'json'], capture_output=True,
                       text=True, timeout=60)
    if r.returncode:
        return None
    s = json.loads(r.stdout).get('status', {})
    if s.get('succeeded'):
        return 'complete'
    if s.get('failed'):
        return 'failed'
    return 'active'


def job_yaml(name, run, shard, gpus=1):
    env = {'COMMIT': 'hybrid-area', 'MANDELBROT_THREADS': str(22 * gpus), 'SHARDED': run['command'],
           'RUN_NAME': run['name'], 'SHARDS': str(run['shards']), 'RUN_TIMEOUT': str(run.get('timeout', 7200)),
           'HOME': '/tmp'}
    if shard is not None:
        env['SHARD'] = str(shard)
    if run.get('build64'):
        env['BUILD64'] = '1'
    envs = ''.join(f'            - name: {k}\n              value: {json.dumps(v)}\n' for k, v in env.items())
    deadline = int(run.get('timeout', 7200)) + 3600
    return f'''apiVersion: batch/v1
kind: Job
metadata:
  name: {name}
  namespace: {NAMESPACE}
  labels:
    kueue.x-k8s.io/queue-name: {QUEUE}
spec:
  activeDeadlineSeconds: {deadline}
  backoffLimit: 0
  ttlSecondsAfterFinished: 604800
  template:
    spec:
      restartPolicy: Never
      securityContext:
        runAsNonRoot: true
        runAsUser: 65532
        runAsGroup: 65532
        fsGroup: 65532
        fsGroupChangePolicy: OnRootMismatch
        seccompProfile:
          type: RuntimeDefault
      volumes:
        - name: data
          persistentVolumeClaim:
            claimName: irving-mandelbrot
        - name: work
          emptyDir:
            sizeLimit: 1Gi
      initContainers:
        - name: fetch
          image: curlimages/curl:8.16.0
          command:
            - sh
            - -c
            - |
              set -e
              mkdir -p /data/bin
              [ -x /data/bin/micromamba ] || {{
                curl -fsSL -o /data/bin/micromamba \\
                  https://github.com/mamba-org/micromamba-releases/releases/download/2.9.0-0/micromamba-linux-64
                chmod +x /data/bin/micromamba; }}
              curl -fsSL -H "Accept: application/vnd.github.raw" -o /work/bench.sh \\
                "https://api.github.com/repos/{REPO}/contents/cluster/bench.sh?ref=hybrid-area"
          volumeMounts:
            - {{ name: data, mountPath: /data }}
            - {{ name: work, mountPath: /work }}
          securityContext:
            allowPrivilegeEscalation: false
            capabilities:
              drop: [ALL]
          resources:
            requests: {{ cpu: "1", memory: 1Gi }}
      containers:
        - name: main
          image: nvcr.io/nvidia/cuda:12.8.1-devel-ubuntu24.04
          env:
{envs}          command: ["bash", "/work/bench.sh"]
          volumeMounts:
            - {{ name: data, mountPath: /data }}
            - {{ name: work, mountPath: /work }}
          securityContext:
            allowPrivilegeEscalation: false
            capabilities:
              drop: [ALL]
          resources:
            requests: {{ cpu: "{22 * gpus}", memory: {64 * gpus}Gi }}
            limits: {{ memory: {64 * gpus}Gi, nvidia.com/gpu: {gpus} }}
'''


def main():
    # Transient kubectl failures (an expired token, API hiccups) should not end a long run: wait and retry
    while True:
        try:
            return loop()
        except (RuntimeError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f'{time.strftime("%H:%M:%S")} kubectl failed ({e}); retrying in 5 minutes '
                  f'(refresh the token with: kubectl auth whoami)', flush=True)
            time.sleep(300)


def loop():
    runs = json.load(open(sys.argv[1]))
    todo = [(run, s) for run in runs for s in range(run['shards'])]
    shard_name = lambda run, s: f'irving-mandelbrot-{run["name"]}-s{s}'
    merge_name = lambda run: f'irving-mandelbrot-{run["name"]}-merge'
    attempts = {}
    merged = set()
    while True:
        states = {(r['name'], s): job_status(shard_name(r, s)) for r, s in todo}
        active = sum(v == 'active' for v in states.values())
        run_active = {r['name']: sum(states[(r['name'], s)] == 'active' for s in range(r['shards'])) for r in runs}
        # Merge runs whose shards are all complete
        for run in runs:
            if run['name'] in merged:
                continue
            if all(states[(run['name'], s)] == 'complete' for s in range(run['shards'])):
                st = job_status(merge_name(run))
                if st is None:
                    kubectl('apply', '-f', '-', input=job_yaml(merge_name(run), run, None))
                    print(f'{time.strftime("%H:%M:%S")} submitted merge of {run["name"]}', flush=True)
                elif st != 'active':
                    merged.add(run['name'])
                    print(f'{time.strftime("%H:%M:%S")} merge of {run["name"]}: {st}; '
                          f'kubectl -n {NAMESPACE} logs job/{merge_name(run)} -c main', flush=True)
        if len(merged) == len(runs):
            return
        # Resubmit failed shards (a few times), and submit new ones while the cluster is idle enough
        queue = None
        for (run, s) in todo:
            key = (run['name'], s)
            if states[key] == 'failed' and attempts.get(key, 0) < 3:
                kubectl('delete', 'job', shard_name(run, s), '--wait=true')
                states[key] = None
            if states[key] is not None:
                continue
            if queue is None:
                queue = queue_state()
            pending, used = queue
            free = run.get('free_gpus', SHORT_FREE_GPUS if run.get('short') else FREE_GPUS)
            if (pending or used > TOTAL_GPUS - free or active >= MAX_JOBS
                    or run_active[run['name']] >= run.get('max_active', MAX_JOBS)):
                continue
            kubectl('apply', '-f', '-', input=job_yaml(shard_name(run, s), run, s))
            attempts[key] = attempts.get(key, 0) + 1
            states[key] = 'active'
            active += 1
            run_active[run['name']] += 1
            print(f'{time.strftime("%H:%M:%S")} submitted {shard_name(run, s)} (queue: {pending} pending, '
                  f'{used} GPUs in use)', flush=True)
            time.sleep(20)  # Let the queue see it before the next check
            queue = None
        time.sleep(60)


if __name__ == '__main__':
    main()

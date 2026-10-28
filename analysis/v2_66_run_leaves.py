"""v2.66: run the leaves of a dumped case split (v2_55_pipeline.py LEAF_DUMP) as separate processes, J at a time,
each with its own certificate file DIR/leaf_K.jsonl and output DIR/leaf_K.txt; a leaf whose certificate already has
its stratum record and all its components is skipped (resumable).
Usage: python3 analysis/v2_66_run_leaves.py QUIVER BLOCK DUMP DIR FIRST LAST [J]"""
import sys, os, json, subprocess, time
Q, B, dump, out, a, b = sys.argv[1:7]; J = int(sys.argv[7]) if len(sys.argv) > 7 else 1
os.makedirs(out, exist_ok=True)
def done(k):
    f = os.path.join(out, f'leaf_{k}.jsonl')
    if not os.path.exists(f): return False
    recs = [json.loads(l) for l in open(f)]
    st = [r for r in recs if r['kind'] == 'stratum']
    return bool(st) and sum(1 for r in recs if r['kind'] != 'stratum') == st[0].get('ncomp', 0)
todo = [k for k in range(int(a), int(b) + 1) if not done(k)]
run = {}
while todo or run:
    while todo and len(run) < J:
        k = todo.pop(0); c = os.path.join(out, f'leaf_{k}.jsonl')
        if os.path.exists(c): os.remove(c)
        env = dict(os.environ, PYTHONHASHSEED='0', BRANCH='1', FOLD_INTEGRAL='1', LEAF_LOAD=dump, LEAVES=str(k), CERT=c)
        run[k] = subprocess.Popen([sys.executable, '-u', 'analysis/v2_55_pipeline.py', Q, B, os.path.join(out, 'scr')],
                                  stdout=open(os.path.join(out, f'leaf_{k}.txt'), 'w'), stderr=subprocess.STDOUT, env=env)
    time.sleep(2)
    for k in list(run):
        if run[k].poll() is not None:
            print(f"leaf {k} rc {run[k].returncode} complete {done(k)}", flush=True); del run[k]

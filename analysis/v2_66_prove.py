"""v2.66: the whole per-quiver proof procedure for Conjecture 3ᶜ (finiteness, uniform in q, every characteristic), run
end to end on one quiver.

Steps (each the script of the round that introduced it):
 1. integral structure: Smith form of Δ, saturated K-basis, det [Kmat | e_C], census group and weightings over ℤ
    (v2_58_integral.py);
 2. vertex blocks: elimination over ℤ of the v2.52 system (P_K conditions and one Δ-condition per parallel class),
    unit leading coefficients, radical certificate, and the point ideal (v2_58_vertex_integral.py, mode int);
 3. char 0 pipeline on every vertex block, with certificates (v2_55_pipeline.py; on failure the integral case split
    BRANCH=1 FOLD_INTEGRAL=1), and the tally (v2_57_tally.py);
 4. specialisation to every prime (v2_58_goodprimes.py, --always 2);
 5. the Bézout bound (v2_58_bound.py, plus the stabiliser factor).
Writes OUT_DIR/summary.json and the logs of every step.
Usage: python3 analysis/v2_66_prove.py QUIVER OUT_DIR"""
import os, sys, json, re, subprocess, collections, math, lzma
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import pcensus as pc
Q, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
SCR = os.path.join(out, 'scratch'); os.makedirs(SCR, exist_ok=True)
def run(args, log, env=None, timeout=20000):
    with open(os.path.join(out, log), 'w') as fh:
        r = subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT, env=dict(os.environ, **(env or {})), timeout=timeout)
    return r.returncode, open(os.path.join(out, log)).read()
summ = {'quiver': Q}
nv, arrows = pc.QUIVERS[Q]
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
# 1
rc, txt = run([sys.executable, 'analysis/v2_58_integral.py', Q], 'integral.txt')
summ['integral'] = dict(snf1='elementary divisors [1]' in txt, saturated='gcd of maximal minors 1' in txt,
                        detB=re.search(r"det B = (-?\d+)", txt).group(1) if 'det B' in txt else None,
                        group_identities=txt.count(': True') >= 6)
# 2
arrow_sel = ','.join(str(c[0]) for c in cls.values())
rc, txt = run([sys.executable, 'analysis/v2_58_vertex_integral.py', Q, os.path.join(SCR, 'vb'), 'int', '--arrow', arrow_sel], 'vertex_integral.txt')
lcs = re.findall(r"^LC (\S+)", txt, re.M); rads = re.findall(r"^RAD \d+ k (\d)", txt, re.M)
m = re.search(r"QDIM (-?\d+) VDIM (\d+)", txt); pts = re.search(r"POINTS (\d+)", txt)
rv = re.search(r"RVDIM (\d+)", txt); rlc = re.search(r"RLC (.*)", txt)
summ['vertex'] = dict(arrows=arrow_sel, zero_dim=bool(m and m.group(1) == '0'), points=int(pts.group(1)) if pts else None,
                      unit_LC=all(x == '1' for x in lcs), radical=all(int(k) > 0 for k in rads) and bool(rads),
                      point_ideal_unit_LC=bool(rlc) and all(x == '1' for x in rlc.group(1).split()),
                      point_ideal_vdim=int(rv.group(1)) if rv else None)
if not summ['vertex']['zero_dim']:
    json.dump(summ, open(os.path.join(out, 'summary.json'), 'w'), indent=1, ensure_ascii=False); print(summ); sys.exit(1)
# the points (vertex blocks), from the char 0 decomposition printed by mode p=0
rc, txt0 = run([sys.executable, 'analysis/v2_58_vertex_integral.py', Q, os.path.join(SCR, 'vb0'), 'p=0', '--arrow', arrow_sel], 'vertex_char0.txt')
blocks = []
for bl in txt0.split("POINTS")[1].split("--")[:-1]:
    vals = {}
    for mm in re.finditer(r"(k\d\d)([-+]\d+)?\s*[,\n]", bl.split("\n", 1)[1] + "\n"):
        vals[mm.group(1)] = -int(mm.group(2)) if mm.group(2) else 0
    assert len(vals) == nv * nv, vals
    blocks.append(','.join(''.join('m' if vals[f'k{i}{j}'] == -1 else str(vals[f'k{i}{j}']) for j in range(nv)) for i in range(nv)))
open(os.path.join(out, 'blocks.txt'), 'w').write("\n".join(blocks) + "\n")
open(os.path.join(out, 'blocks_census_format.txt'), 'w').write('vertex block |\n' + ''.join(
    ' '.join(r.replace('m', '-1') for r in b.split(',')) + ' | 0 0 0\n' for b in blocks))
# 3
cert = os.path.join(out, 'certs.jsonl')
if os.path.exists(cert): os.remove(cert)
plog = open(os.path.join(out, 'pipeline.txt'), 'w'); reruns = []
for b in blocks:
    for mode in ({}, {'BRANCH': '1', 'FOLD_INTEGRAL': '1'}):
        tmpc = cert + '.tmp'
        if os.path.exists(tmpc): os.remove(tmpc)
        try:
            r = subprocess.run([sys.executable, '-u', 'analysis/v2_55_pipeline.py', Q, b, os.path.join(SCR, 'p_' + b)], cwd=ROOT,
                               env=dict(os.environ, CERT=tmpc, **mode), capture_output=True, text=True, timeout=6000)
            rc_, o = r.returncode, r.stdout + r.stderr
        except subprocess.TimeoutExpired: rc_, o = -1, 'TIMEOUT'
        if rc_ == 0:
            plog.write("\n".join(l for l in o.splitlines() if not l.startswith('leaf ')) + "\n"); plog.flush()
            with open(cert, 'a') as fh: fh.write(open(tmpc).read())
            if mode: reruns.append(b)
            break
    else:
        plog.write(f"### {b}: FAILED in both modes\n"); plog.flush()
plog.close()
rc, txt = run([sys.executable, 'analysis/v2_57_tally.py', Q, os.path.join(out, 'blocks_census_format.txt'), os.path.join(out, 'pipeline.txt')], 'tally.txt')
last = txt.strip().splitlines()[-1]
summ['pipeline'] = dict(tally=last, case_split_blocks=reruns)
# pivots / constants
piv = collections.Counter(); bc = set()
for l in open(cert):
    r = json.loads(l)
    if r['kind'] == 'stratum':
        for _, c in r.get('pivots', []): piv[c] += 1
        bc |= set(map(tuple, r.get('branch_consts', [])))
summ['pivots'] = dict(piv); summ['branch_consts'] = sorted(bc)
# 4
rc, txt = run([sys.executable, 'analysis/v2_58_goodprimes.py', cert, os.path.join(SCR, 'gp'), '--timeout', '2400', '--always', '2'], 'goodprimes.txt')
cands = collections.Counter(); claims = 0; strata = 0
for l in txt.splitlines():
    mm = re.match(r"^(\S+) leaf (\d+): \d+ components, (\d+) specialisation claims; candidate primes \[(.*?)\]", l)
    if mm:
        strata += 1; claims += int(mm.group(3))
        for p in mm.group(4).split(','):
            if p.strip(): cands[int(p)] += 1
exc = re.search(r"EXCEPTIONAL (.*)", txt)
summ['goodprimes'] = dict(strata=strata, claims=claims, candidates=dict(cands), exceptional=exc.group(1) if exc else 'unresolved',
                          cross_check_unresolved=txt.count('cross-check unresolved'))
# 5
rc, txt = run([sys.executable, 'analysis/v2_58_bound.py', cert, Q], 'bound.txt')
lb = float(re.search(r"10\^(-?[\d.]+)", txt.strip().splitlines()[-1]).group(1))
maxlen = max(len(b[1]) for b in pc.setup(Q, 3)['basis'] if b[0] == 'p')
Ncoord = sum(len(c) * (len(c) - 1) + 1 for c in cls.values() if len(c) > 1) + 2 * len(cls)
summ['bound_log10'] = round(lb + Ncoord * math.log10(2 * maxlen), 1)
json.dump(summ, open(os.path.join(out, 'summary.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps(summ, ensure_ascii=False))

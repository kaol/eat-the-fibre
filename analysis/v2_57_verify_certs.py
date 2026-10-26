"""v2.57: check the orbit-dimension certificates written by v2_55_pipeline.py (env CERT=file.jsonl).

Per stratum record (one per vertex block, or per leaf in BRANCH mode) and each of its component records:
 (a) P is rebuilt from fold_ideal.stratum(QUIVER, V) and the recorded substitutions; 𝔥_Z from the recorded Der_Δ basis,
     class scalings and W_Z; and every recorded minor is recomputed as the determinant of the recorded rows (elements
     of 𝔥_Z) and columns (entries of P) of M = ([h, P])_h, and compared with the recorded polynomial;
 (b) Singular: the presolved equations of the stratum (and, in BRANCH mode, nothing more: saturation only removes
     components) reduce to 0 modulo the component's prime, i.e. Z lies in the stratum;
 (c) Singular, with std and dim only (no primary decomposition): for every r < dim Z, dim(prime + recorded minors) ≤ r.
     Since the minors vanish on Z_r, this is dim Z_r ≤ r.
With --coverage, also (d): every generator of the intersection of the primes of a stratum record lies in the radical of
its ideal (saturated by the inequations in BRANCH mode), so the primes cover it (Rabinowitsch, per generator; timeout).
Usage: python3 analysis/v2_57_verify_certs.py CERT.jsonl[.xz] [--coverage] [--timeout S]"""
import os, sys, json, subprocess, argparse, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import fold_ideal as fi
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('--coverage', action='store_true')
ap.add_argument('--timeout', type=int, default=3600); ap.add_argument('--tmp', default='/tmp')
A = ap.parse_args()
import lzma
recs = [json.loads(l) for l in (lzma.open(A.cert, 'rt') if A.cert.endswith('.xz') else open(A.cert))]
groups = collections.OrderedDict()
for r in recs:
    key = (r['quiver'], r['block'], r['leaf'])
    g = groups.setdefault(key, {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
SS = {}
def stratum(q, b):
    if (q, b) not in SS:
        SS[(q, b)] = fi.stratum(q, [[-1 if ch == 'm' else int(ch) for ch in row] for row in b.split(',')])
    return SS[(q, b)]
ss = lambda f: str(f).replace('**', '^')
def singular(src, tag):
    fn = os.path.join(A.tmp, f'vc_{os.getpid()}_{tag}.sing'); open(fn, 'w').write(src)
    try: res = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=A.timeout)
    except subprocess.TimeoutExpired: return None
    if "? error" in res.stdout + res.stderr: return None
    return res.stdout
tot = collections.Counter()
for (q, b, leaf), g in groups.items():
    st = g['stratum']; comps = g['comps']
    if st is None: print(f"{b} leaf {leaf}: no stratum record"); tot['missing'] += 1; continue
    S = stratum(q, b); R = S['R']; d = S['ctx']['d']
    gl = {str(x): x for x in R.gens}
    parse = lambda s_: R(sp.sympify(s_, locals={k: sp.Symbol(k) for k in gl}))
    vs = st['vs']
    if len(comps) != st['ncomp']: print(f"{b} leaf {leaf}: {len(comps)} component records, {st['ncomp']} expected"); tot['missing'] += 1
    if not vs:
        tot['point strata'] += 1; continue
    # (a) minors
    Pr = [row[:] for row in S['P']]
    for x, v in st['subs']:
        xv, vv = gl[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    DerD = [sp.Matrix(d, d, [sp.Rational(x) for x in D]) for D in st['DerD']]
    Ec = [sp.diag(*e) for e in st['Ec']]
    KD = R.to_domain(); ok_a = True
    for c in comps:
        if not c['cuts']: continue
        W = [[sp.Rational(x) for x in w] for w in c['W']]
        h = DerD + [sum((w[k] * Ec[k] for k in range(len(Ec))), sp.zeros(d, d)) for w in W]
        Mcache = {}
        def Mrow(i):
            if i not in Mcache:
                D = h[i]
                Dr = [[R(sp.QQ(int(sp.fraction(D[a, k])[0]), int(sp.fraction(D[a, k])[1]))) if D[a, k] else R(0) for k in range(d)] for a in range(d)]
                Mcache[i] = [sum((Dr[a][k] * Pr[k][bb] for k in range(d) if Dr[a][k]), R(0)) -
                             sum((Pr[a][k] * Dr[k][bb] for k in range(d) if Dr[k][bb]), R(0)) for a in range(d) for bb in range(d)]
            return Mcache[i]
        for cut in c['cuts']:
            for mn in cut['minors']:
                sub = [[Mrow(i)[j] for j in mn['cols']] for i in mn['rows']]
                det = DomainMatrix(sub, (len(sub), len(sub)), KD).det()
                if len(mn['rows']) != cut['r'] + 1 or str(det) != mn['det']:
                    ok_a = False; print(f"{b} leaf {leaf} comp {c['comp']} r {cut['r']}: minor mismatch")
    # (b), (c) in one Singular run
    lines = [f"ring R = 0,({','.join(vs)}),dp;", f"ideal I = {','.join(ss(e) for e in st['rem']) or '0'};", "int bad = 0;"]
    for c in comps:
        if c['dim'] <= 0 and not c['gens']: continue
        lines.append(f"ideal G = {','.join(ss(x) for x in c['gens']) or '0'}; ideal SG = std(G);")
        lines.append(f'if (size(reduce(I, SG)) != 0) {{ print("NOTIN {c["comp"]}"); }}')
        lines.append(f'if (dim(SG) != {c["dim"]}) {{ print("DIMWRONG {c["comp"]}"); }}')
        for cut in c['cuts']:
            F = [ss(m['det']) for m in cut['minors']]
            lines.append(f"ideal J = SG{',' if F else ''}{','.join(F)}; "
                         f'if (dim(std(J)) > {cut["r"]}) {{ print("CUTFAIL {c["comp"]} {cut["r"]}"); }} else {{ print("CUTOK {c["comp"]} {cut["r"]}"); }} kill J;')
        lines.append("kill G; kill SG;")
    if A.coverage and comps:
        lines.append('LIB "elim.lib";')
        lines += [f"I = sat(I, {ss(f)});" for f in st['nz']]
        lines.append("ideal K = " + ("intersect(" + ",".join(f"ideal({','.join(ss(x) for x in c['gens']) or '0'})" for c in comps) + ")" if len(comps) > 1
                                     else f"ideal({','.join(ss(x) for x in comps[0]['gens']) or '0'})") + ";")
        lines.append(f"ring R2 = 0,({','.join(vs)},TT),dp; ideal I2 = imap(R, I); ideal K2 = imap(R, K); int cv = 1; int ii;")
        lines.append('for (ii = 1; ii <= size(K2); ii++) { ideal T = std(I2 + ideal(1 - TT * K2[ii])); '
                     'if (!(size(T) == 1 && deg(T[1]) == 0)) { cv = 0; } kill T; } print("COVER " + string(cv));')
    lines.append('print("DONE"); quit;')
    out = singular("\n".join(lines), f"{b.replace(',', '')}_{leaf}")
    if out is None or "DONE" not in out:
        print(f"{b} leaf {leaf}: Singular timeout/error"); tot['singular unresolved'] += 1; continue
    bad = [l for l in out.splitlines() if l.startswith(("NOTIN", "DIMWRONG", "CUTFAIL"))]
    ncut = sum(1 for l in out.splitlines() if l.startswith("CUTOK"))
    cov = [l for l in out.splitlines() if l.startswith("COVER")]
    tot['strata'] += 1; tot['components'] += len(comps); tot['cuts ok'] += ncut; tot['bad'] += len(bad) + (0 if ok_a else 1)
    if A.coverage: tot['coverage ' + (cov[0].split()[1] if cov else '?')] += 1
    print(f"{b} leaf {leaf}: {len(comps)} components, {ncut} cuts verified, minors {'ok' if ok_a else 'MISMATCH'}"
          + (f", {bad}" if bad else "") + (f", {cov[0]}" if cov else ""), flush=True)
print("TOTAL", dict(tot))

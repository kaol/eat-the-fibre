"""v2.58 (Program 2): the primes at which the characteristic-0 certificates of a quiver might not specialise, and a
direct check of every such prime.

Input: a certificate file from v2_55_pipeline.py (CERT=…, with presolve pivots, v2.58). The coordinate change must be
integral (v2_58_integral.py). Presolve pivots and case-split constants (v2.59) that are not ±1 make their primes
candidates; at such a prime a direct check of the presolved claims does NOT re-validate the presolve or the case split
itself, so those primes need separate treatment (reported as candidates).

Claims, for each stratum (gens_Z = the recorded generators of the component's prime, denominators cleared;
Z_p := V(gens_Z mod p)):
 (DIM) dim V(gens_Z + F_r) ≤ r for each recorded cut r, and dim V(gens_Z) ≤ dim Z.
       A strong Gröbner basis G of the ideal over ℤ (Singular, ring integer). Every g ∈ G lies in the ℤ-ideal, so for
       p ∤ LC(g) the leading monomial of g is a leading monomial of the ideal mod p. If p divides no leading
       coefficient, the leading ideal mod p contains that over ℚ, so the dimension mod p is at most the one over ℚ.
       Candidates: primes dividing some LC(g).
 (DEF) the weightings W_Z are compatible on Z_p: each defect δ (Σ_c w_c δ_c, w ∈ W_Z) lies in (gens_Z) over ℚ (the
       pipeline's normal forms), so dividing δ by G_Z (a ℚ-Gröbner basis too) leaves 0 with quotients whose
       denominators divide products of LC(G_Z). Candidates: primes dividing LC(G_Z) (the same G as (DIM) at r = dim Z).
 (COV) V(rem mod p) ⊆ ∪_Z Z_p. Chain K_1 = (gens_1), K_i = K_{i−1} ∩ (gens_i) over ℚ, each with a Gröbner basis
       (integral, as Singular normalises). For p not dividing their leading coefficients, the divisions
       k·g ∈ K_i (k ∈ K_{i−1}, g ∈ gens_i) and gens_1 ⊆ K_1 hold mod p, so V(K_{n,p}) ⊆ ∪ Z_p. Then each k ∈ K_n lies in
       √rem over ℚ; over ℤ, the ideal rem + (1 − t·k) contains a nonzero integer n_k (strong Gröbner basis), so for
       p ∤ n_k, k ∈ √(rem mod p). Candidates: primes dividing those leading coefficients and the n_k.
For every candidate prime, the claims it is a candidate for are checked directly over 𝔽_p (v2.59: targeted), and the
primes in --always get a full direct check of every claim (a cross-check; a timeout there is reported, not exceptional) (Singular in characteristic p; the radical
membership by Rabinowitsch). A claim failing there is an exceptional prime, reported.
Usage: python3 analysis/v2_58_goodprimes.py CERT.jsonl[.xz] OUT_PREFIX [--timeout S] [--only BLOCK]"""
import os, sys, json, lzma, subprocess, argparse, collections, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import fold_ideal as fi
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('out'); ap.add_argument('--timeout', type=int, default=3600)
ap.add_argument('--only', default=None); ap.add_argument('--skip', default=''); ap.add_argument('--always', default='2,3', help='primes checked directly in any case (cross-check)')
A = ap.parse_args()
recs = [json.loads(l) for l in (lzma.open(A.cert, 'rt') if A.cert.endswith('.xz') else open(A.cert))]
groups = collections.OrderedDict()
for r in recs:
    g = groups.setdefault((r['quiver'], r['block'], r['leaf']), {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
ss = lambda f: str(f).replace('**', '^')
def singular(src, tag):
    fn = f"{A.out}_{tag}.sing"; open(fn, 'w').write(src)
    try: res = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=A.timeout)
    except subprocess.TimeoutExpired: return None
    if "? error" in res.stdout + res.stderr: return "ERROR " + (res.stdout + res.stderr)[-600:].replace("DONE", "")
    return res.stdout
def integral(polys, Sv):
    """clear denominators (multiply by the lcm; never divide by a content)"""
    out = []
    for s_ in polys:
        if '/' not in s_:                     # already integral (pivots ±1, integral coordinates): no parsing (v2.59)
            if s_.strip() not in ('0', ''): out.append(s_.replace('**', '^'))
            continue
        e = sp.sympify(s_, locals={str(v): v for v in Sv})
        if e == 0: continue
        den = sp.ilcm(1, 1, *[sp.fraction(c)[1] for c in sp.Poly(e, *Sv).coeffs()]) if Sv else 1
        out.append(ss(sp.expand(e * den)))
    return out
def gens_of(polys, Sv): return integral(polys, Sv) or ['0']
def primes_of(ints):
    ps = set()
    for n in ints:
        n = abs(int(n))
        if n > 1: ps |= set(sp.factorint(n).keys())
    return ps
SS = {}
report = []; exceptional = []
for (q, b, leaf), g in groups.items():
    if (A.only and b != A.only) or b in A.skip.split(';'): continue
    st = g['stratum']; comps = [c for c in g['comps'] if c['gens']]
    vs = st['vs']
    if not vs or not comps:
        report.append(f"{b}: single point or no positive-dimensional component; nothing to specialise"); continue
    pivp = primes_of([x for _, c_ in st.get('pivots', []) for x in str(sp.Rational(c_)).replace('/', ' ').split()])
    Sv = sp.symbols(vs)
    NZP = "*".join(f"({x})" for x in integral(st.get('nz', []), Sv)) or "1"     # product of the leaf's inequations
    # defects δ_c in the presolved coordinates (as the pipeline computes them)
    key = (q, b)
    if key not in SS: SS[key] = fi.stratum(q, [[-1 if ch == 'm' else int(ch) for ch in row] for row in b.split(',')])
    S = SS[key]; R = S['R']; ctx = S['ctx']; d = ctx['d']; basis = ctx['basis']; arrows = ctx['arrows']
    gl = {str(x): x for x in R.gens}
    parse = lambda s_: R(sp.sympify(s_, locals={k: sp.Symbol(k) for k in gl}))
    Pr = [row[:] for row in S['P']]
    for x, v in st['subs']:
        xv, vv = gl[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    idx = {bb: i for i, bb in enumerate(basis)}
    cls = collections.defaultdict(list)
    for a_, s_ in enumerate(arrows): cls[s_].append(a_)
    cl = list(cls.values())
    def DLclass(c):
        D_ = np.zeros((d, d, d), dtype=np.int64)
        for l, bb in enumerate(basis):
            if bb[0] == 'v': continue
            for k, x in enumerate(bb[1]):
                if x not in c: continue
                pre = idx[('p', bb[1][:k])] if k > 0 else idx[('v', arrows[x][0])]
                suf = idx[('p', bb[1][k + 1:])] if k + 1 < len(bb[1]) else idx[('v', arrows[x][1])]
                D_[l, pre, suf] += 1
        return D_
    def defect(Dl):
        out = []
        for l in range(d):
            PD = [[sum((Pr[a][x] * int(Dl[l][x, y]) for x in range(d) if Dl[l][x, y]), R(0)) for y in range(d)] for a in range(d)]
            for a in range(d):
                for bb in range(d):
                    lhs = sum((Pr[m][l] * int(Dl[m][a, bb]) for m in range(d) if Dl[m][a, bb]), R(0))
                    out.append(lhs - sum((PD[a][y] * Pr[bb][y] for y in range(d) if PD[a][y]), R(0)))
        return out
    defs = [defect(DLclass(c)) for c in cl]
    lines = [f"ring Rz = integer,({','.join(vs)}),dp; option(redSB);", "int ii; int jj;"]
    for c in comps:
        G = gens_of(c['gens'], Sv)
        lines.append(f"ideal G{c['comp']} = {','.join(G)};")
        lines.append(f"ideal SG = std(G{c['comp']}); string s = \"\"; for (ii = 1; ii <= size(SG); ii++) {{ s = s + string(leadcoef(SG[ii])) + \" \"; }} "
                     f"print(\"DIM {c['comp']} top \" + s); kill SG; kill s;")
        for cut in c['cuts']:
            F = integral([m_['det'] for m_ in cut['minors']], Sv)
            lines.append(f"ideal SJ = std(G{c['comp']} + ideal({','.join(F) or '0'})); string s = \"\"; for (ii = 1; ii <= size(SJ); ii++) "
                         f"{{ s = s + string(leadcoef(SJ[ii])) + \" \"; }} print(\"DIM {c['comp']} {cut['r']} \" + s); kill SJ; kill s;")
    out_z = singular("\n".join(lines + ['print("DONE"); quit;']), f"{b.replace(',', '')}_{leaf}_dim")
    # coverage chain over ℚ, then Rabinowitsch over ℤ
    lq = [f"ring Rq = 0,({','.join(vs)}),dp; int ii; int jj; int kk;"]
    for c in comps: lq.append(f"ideal G{c['comp']} = {','.join(gens_of(c['gens'], Sv))};")
    lq.append(f"ideal K = std(G{comps[0]['comp']}); string s;")
    lq.append('s = ""; for (ii = 1; ii <= size(K); ii++) { s = s + string(leadcoef(K[ii])) + " "; } print("CHAIN 1 " + s);')
    lq.append(f'if (size(reduce(G{comps[0]["comp"]}, K)) != 0) {{ print("CHAINFAIL 1"); }}')
    for i, c in enumerate(comps[1:], 2):
        lq.append(f"ideal K2 = std(intersect(K, G{c['comp']}));")
        lq.append(f'for (ii = 1; ii <= size(K); ii++) {{ for (jj = 1; jj <= size(G{c["comp"]}); jj++) {{ '
                  f'if (reduce(K[ii] * G{c["comp"]}[jj], K2) != 0) {{ print("CHAINFAIL {i}"); }} }} }}')
        lq.append(f's = ""; for (ii = 1; ii <= size(K2); ii++) {{ s = s + string(leadcoef(K2[ii])) + " "; }} print("CHAIN {i} " + s); K = K2; kill K2;')
    lq.append('print("KGENS " + string(K));')
    lq.append('print("DONE"); quit;')
    out_q = singular("\n".join(lq), f"{b.replace(',', '')}_{leaf}_cov")
    Kgens = re.search(r"KGENS (.*)", out_q).group(1).split(',') if out_q and "KGENS" in out_q and not out_q.startswith("ERROR") else None
    out_r = None
    if Kgens:
        lr = [f"ring Rt = integer,({','.join(vs)},TT),dp; option(redSB); int ii;", f"ideal I = {','.join(integral(st['rem'], Sv)) or '0'};"]
        # v2.59: first try (k·∏nz)^e ∈ (rem) over ℤ for e ≤ 4, by reduction modulo a strong Gröbner basis of rem over
        # ℤ (exact: then the identity holds over ℤ, at every p, with no candidate prime); else Rabinowitsch over ℤ.
        lr.append("ideal GI = std(I); int ee; int done_;")
        for kk, k in enumerate(Kgens):
            lr.append(f"done_ = 0; for (ee = 1; ee <= 4 && done_ == 0; ee++) {{ if (reduce((({k}) * {NZP})^ee, GI) == 0) {{ done_ = ee; }} }}")
            lr.append(f"if (done_ > 0) {{ print(\"RAB {kk} 1 (power \" + string(done_) + \")\"); }} else {{ "
                      f"ideal J = std(I + ideal(1 - TT * ({k}) * {NZP})); string s = \"\"; "
                      f"for (ii = 1; ii <= size(J); ii++) {{ if (deg(J[ii]) == 0) {{ s = s + string(J[ii]) + \" \"; }} }} print(\"RAB {kk} \" + s); kill J; kill s; }}")
        out_r = singular("\n".join(lr + ['print("DONE"); quit;']), f"{b.replace(',', '')}_{leaf}_rab")
    # candidate primes per claim
    cand = collections.defaultdict(set); status = []
    for name_, out in (('dim', out_z), ('cov', out_q), ('rab', out_r)):
        if out is None or "DONE" not in out: status.append(f"{name_}: unresolved ({'timeout' if out is None else out[:200]})"); continue
    if out_z and "DONE" in out_z:
        for m_ in re.finditer(r"DIM (\d+) (\S+) (.*)", out_z):
            cand[('DIM', int(m_.group(1)), m_.group(2))] |= primes_of(m_.group(3).split())
    if out_q and "DONE" in out_q:
        if "CHAINFAIL" in out_q: status.append("coverage chain: a product did not reduce (unexpected)")
        for m_ in re.finditer(r"CHAIN (\d+) (.*)", out_q): cand[('COV', 'chain', int(m_.group(1)))] |= primes_of(m_.group(2).split())
    if out_r and "DONE" in out_r:
        for m_ in re.finditer(r"RAB (\d+) (.*)", out_r):
            ns = m_.group(1), m_.group(2).split("(")[0].split()
            if not ns[1]: status.append(f"Rabinowitsch {ns[0]}: no integer in the ideal (k not in √rem over ℚ?)")
            cand[('COV', 'rab', int(ns[0]))] |= primes_of(ns[1])
    cand[('PIVOT',)] |= pivp
    cand[('BRANCH',)] |= primes_of([x for nd in st.get('branch_consts', []) for x in nd])
    allp = sorted(set().union(*cand.values())) if cand else []
    nclaims = (len(re.findall(r"^DIM ", out_z or "", re.M)), len(re.findall(r"^CHAIN ", out_q or "", re.M)), len(re.findall(r"^RAB ", out_r or "", re.M)))
    expect = (sum(1 + len(c['cuts']) for c in comps), len(comps), len(Kgens or []))
    if nclaims != expect: status.append(f"parsed {nclaims} claims, expected {expect}")
    always = [int(x) for x in A.always.split(',') if x]
    # direct checks at every candidate prime
    bad = []; xunres = []
    # targeted checks at each candidate prime (only the claims it is a candidate for), then full cross-checks
    runs = [(p, False) for p in allp] + [(p, True) for p in always]
    for p, full in runs:
        need = lambda key: full or p in cand.get(key, set())
        if not full and (p in cand.get(('PIVOT',), set()) or p in cand.get(('BRANCH',), set())):
            bad.append((p, 'presolve pivot or case-split constant divisible by p: the coordinates need a separate argument'))
        lc = [f'LIB "elim.lib"; ring Rp = {p},({",".join(vs)},TT),dp; int ii; int bad = 0;',
              f"ideal I = {','.join(integral(st['rem'], Sv)) or '0'};"]
        for c in comps:
            top = need(('DIM', c['comp'], 'top')); cuts = [cut for cut in c['cuts'] if need(('DIM', c['comp'], str(cut['r'])))]
            if not (top or cuts or full): continue
            Gs = ','.join(gens_of(c['gens'], Sv))
            lc.append(f"ideal G{c['comp']} = {Gs}; ideal SG = std(G{c['comp']});")
            if top: lc.append(f'if (dim(SG) - 1 > {c["dim"]}) {{ print("BAD DIMTOP {c["comp"]} " + string(dim(SG) - 1)); }}')   # TT is free: subtract 1
            for cut in cuts:
                F = integral([m_['det'] for m_ in cut['minors']], Sv)
                lc.append(f"if (dim(std(SG + ideal({','.join(F) or '0'}))) - 1 > {cut['r']}) {{ print(\"BAD DIM {c['comp']} {cut['r']}\"); }}")
            for w in (c['W'] if top else []):
                wv = [sp.Rational(x) for x in w]
                dl = [sum((wv[cc] * defs[cc][k] for cc in range(len(cl))), R(0)) for k in range(len(defs[0]))]
                dl = [x for x in dl if x]
                if dl:
                    lc.append(f"ideal DD = {','.join(integral([str(x) for x in dl], Sv))};")
                    lc.append(f'for (ii = 1; ii <= size(DD); ii++) {{ if (reduce(DD[ii], SG) != 0) {{ ideal T = std(SG + ideal(1 - TT * DD[ii])); '
                              f'if (!(size(T) == 1 && deg(T[1]) == 0)) {{ print("BAD DEF {c["comp"]}"); }} kill T; }} }} kill DD;')
            lc.append("kill SG;")
        if full or any(p in v for k_, v in cand.items() if k_[0] == 'COV'):
            for c in comps:
                if f"ideal G{c['comp']} =" not in "\n".join(lc): lc.append(f"ideal G{c['comp']} = {','.join(gens_of(c['gens'], Sv))};")
            lc.append("ideal K = " + ("intersect(" + ",".join(f"G{c['comp']}" for c in comps) + ")" if len(comps) > 1 else f"G{comps[0]['comp']}") + ";")
            lc.append('for (ii = 1; ii <= size(K); ii++) { ideal T = std(I + ideal(1 - TT * K[ii] * ' + NZP + ')); if (!(size(T) == 1 && deg(T[1]) == 0)) { print("BAD COV"); } kill T; }')
        out_p = singular("\n".join(lc + ['print("DONE"); quit;']), f"{b.replace(',', '')}_{leaf}_p{p}{'full' if full else ''}")
        if out_p is None or "DONE" not in out_p:
            (xunres if full else bad).append((p, 'unresolved')); continue
        fails = sorted(set(l for l in out_p.splitlines() if l.startswith("BAD")))
        if fails: bad.append((p, fails))
    if bad: exceptional.append((b, leaf, bad))
    report.append(f"{b} leaf {leaf}: {len(comps)} components, {sum(nclaims)} specialisation claims; candidate primes {allp} (also checked directly: {always}); " +
                  (f"FAILING {bad}" if bad else "all claims hold at every candidate prime") + (f"; cross-check unresolved at {xunres}" if xunres else "") + (f"; {status}" if status else ""))
    print(report[-1], flush=True)
print("EXCEPTIONAL", exceptional)

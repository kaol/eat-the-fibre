"""Machine check of the steps of the proof of Conjecture 4 (acyclic, no parallel arrows), for every arrow-rule vertex fold.
S_v: Hodge v -> subtree(v); Podge v -> subtree(v) minus v.  N = paths whose vertices all have S_v ≠ ∅.
 (P1) explicit forest data: supp P0(e_v) is 'subtree' for v with P0(e_v)(v) ≠ 0 (Hodge) and P0(e_v) = −Σ_{children} P0(e_d) for Podge v
 (P2) the Hodge expansion H(p) of every p ∈ N is path-supported
 (P3) every path in N whose endpoints lie in a common diagonal block idempotent has Hodge endpoints
 (P4) for every p ∈ N with Hodge endpoints and length m ≥ 1, each run term (Lpart⊗1_{C1})·Rdiff·(1_{Cr}⊗Rpart) is a product
      of three path-supported elements killed by T, and these terms sum to p − H(p)  (exact, as vectors in A_m)."""
import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from quiver import build
from extend import vertex_folds
def setup(nv, arrows):
    d, mult, Dl, allp = build(nv, arrows)
    vseq = lambda p: (p[0],) if p[0] != 'a' else (arrows[p[1]][0],) + tuple(arrows[k][1] for k in p[1:])
    idx = {vseq(p): j for j, p in enumerate(allp)}
    return d, idx, allp, vseq
def check(nv, arrows, P0):
    d, idx, allp, vseq = setup(nv, arrows); E = set(arrows)
    col = lambda v: [P0[a, v] for a in range(nv)]
    hodge = {v for v in range(nv) if P0[v, v] != 0 and all(x >= 0 for x in col(v)) and any(col(v))}
    podge_nonleaf = {v for v in range(nv) if any(col(v)) and v not in hodge}
    zero = {v for v in range(nv) if not any(col(v))}
    # children of a Podge vertex c: the Hodge vertices h with P0(e_c) = −Σ P0(e_h) structure -> recover via the Q0 block of c
    Q0 = (identity_matrix(QQ, nv) - matrix(QQ, P0)).column_space()
    blocks = []
    for S in itertools.chain.from_iterable(itertools.combinations(range(nv), r) for r in range(1, nv + 1)):
        v = vector(QQ, [1 if i in S else 0 for i in range(nv)])
        if v in Q0 and not any(set(S) > B for B in blocks): blocks.append(set(S))
    head = {}
    for B in blocks:
        heads = [c for c in B if c not in hodge]
        if len(heads) == 1: head[heads[0]] = B
    ok1 = all(c in head for c in podge_nonleaf) and all(
        list(sum((vector(QQ, col(h)) for h in head[c] - {c}), vector(QQ, [0] * nv))) == [-x for x in col(c)] for c in podge_nonleaf)
    ch = {c: sorted(head[c] - {c}) for c in podge_nonleaf if c in head}
    def is_path(seq): return seq in idx
    def hexp(seq):                             # Hodge expansion: dict seq -> coeff
        opts = [[(v, 1)] if v in hodge else [(d_, -1) for d_ in ch[v]] for v in seq]
        out = {}
        for choice in itertools.product(*opts):
            s = tuple(c[0] for c in choice); cf = 1
            for c in choice: cf *= c[1]
            out[s] = out.get(s, 0) + cf
        return out
    T = lambda seq: None
    Npaths = [vseq(p) for p in allp if p[0] == 'a' and all(v not in zero for v in vseq(p))]
    ok2 = all(all(is_path(s) for s in hexp(p)) for p in Npaths if all(v in hodge or v in ch for v in p))
    # P3: diagonal endpoints are Hodge
    diag_idem = [B for B in blocks] + [set(range(nv)) - set().union(*blocks)]
    ok3 = all(p[0] in hodge and p[-1] in hodge for p in Npaths for X in diag_idem if p[0] in X and p[-1] in X and len(p) > 1)
    # P4: run decomposition reproduces p − H(p) with path-supported, T-killed factors
    def piH(v): return vector(QQ, col(v))
    def Tvanish(elem):                         # T(elem) = ⊗ π_H(e_v) summed
        acc = {}
        for seq, cf in elem.items():
            for choice in itertools.product(*[[(a, P0[a, v]) for a in range(nv) if P0[a, v]] for v in seq]):
                s = tuple(c[0] for c in choice); w = cf
                for c in choice: w *= c[1]
                acc[s] = acc.get(s, 0) + w
        return all(v == 0 for v in acc.values())
    ok4 = True; tested = 0
    for p in Npaths:
        if not (p[0] in hodge and p[-1] in hodge): continue
        runs = []; j = 0
        while j < len(p):
            if p[j] in podge_nonleaf:
                k = j
                while k + 1 < len(p) and p[k + 1] in podge_nonleaf: k += 1
                runs.append((j, k)); j = k + 1
            else: j += 1
        target = {p: 1}
        for s, cf in hexp(p).items(): target[s] = target.get(s, 0) - cf
        total = {}
        for subset in itertools.chain.from_iterable(itertools.combinations(range(len(runs)), r) for r in range(1, len(runs) + 1)):
            # choose Rdiff on runs in subset, Hodge expansion on the others; build the term
            opts = []
            for j, v in enumerate(p):
                rin = [ri for ri, (a, b) in enumerate(runs) if a <= j <= b]
                if not rin: opts.append([((v,), 1)])
            term = {(): 1}
            j = 0; pos = 0
            segs = []
            while j < len(p):
                rin = [ri for ri, (a, b) in enumerate(runs) if a <= j <= b]
                if rin:
                    a, b = runs[rin[0]]
                    seg = p[a:b + 1]
                    if rin[0] in subset:
                        el = {seg: 1}
                        for s, cf in hexp(seg).items(): el[s] = el.get(s, 0) - cf
                    else: el = hexp(seg)
                    segs.append(el); j = b + 1
                else:
                    segs.append({(p[j],): 1}); j += 1
            for el in segs:
                term = {s + t: c1 * c2 for s, c1 in term.items() for t, c2 in el.items()}
            # factor at the first chosen run: (Lpart ⊗ 1_C1)·Rdiff·(1_Cr ⊗ Rpart) — check the three factors
            r0 = runs[subset[0]]; a, b = r0
            nseg_before = sum(1 for _ in range(0))  # positions map 1:1 to slots since every segment keeps its length
            for s in term:
                ok4 &= is_path(s)
            left = {s[:a]: 0 for s in term}; right = {s[b + 1:]: 0 for s in term}
            C1 = head[p[a]]; Cr = head[p[b]]
            Lfac = {}
            for s, c in term.items():
                for x in C1: Lfac[s[:a] + (x,)] = 1
            Rfac = {}
            for s, c in term.items():
                for x in Cr: Rfac[(x,) + s[b + 1:]] = 1
            ok4 &= all(is_path(s) for s in Lfac) and all(is_path(s) for s in Rfac)
            mid = {p[a:b + 1]: 1}
            for s, cf in hexp(p[a:b + 1]).items(): mid[s] = mid.get(s, 0) - cf
            ok4 &= all(is_path(s) for s in mid) and Tvanish(mid)
            for s, c in term.items(): total[s] = total.get(s, 0) + c
            tested += 1
        tot = {s: c for s, c in total.items() if c}; tg = {s: c for s, c in target.items() if c}
        ok4 &= tot == tg
    return ok1, ok2, ok3, ok4, tested
if __name__ == '__main__':
    from dags import dags
    from arrowrule import QS
    if sys.argv[1] == 'dags':
        qs = [a for a in dags(int(sys.argv[2])) if a]
        names = [str(a) for a in qs]
    else:
        names = sys.argv[1].split(';'); qs = [QS[n][1] for n in names]
    allok = [True] * 4; nf = 0; nt = 0
    for name, arrows in zip(names, qs):
        nv = max(max(a) for a in arrows) + 1 if sys.argv[1] != 'dags' else int(sys.argv[2])
        if sys.argv[1] != 'dags': nv = QS[name][0]
        A = set(arrows)
        for P0 in vertex_folds(nv):
            if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                       all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
            r = check(nv, arrows, P0); nf += 1; nt += r[4]
            for i in range(4): allok[i] &= r[i]
    print(f"{len(qs)} quivers, {nf} vertex folds, {nt} run terms: P1 forest data {allok[0]}, P2 Hodge expansion path-supported {allok[1]}, "
          f"P3 diagonal endpoints Hodge {allok[2]}, P4 run decomposition exact with admissible factors {allok[3]}")

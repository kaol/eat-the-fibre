"""Are canonical extensions and loop-collapse folds all the folds of a looped quiver without parallel arrows?
For every arrow-rule vertex fold, compare the exact truncated fibre dimension with the collapse prediction:
  P0 is a collapse vertex fold iff P0 e_s = e_s (s looped), P0 e_{s'} = −e_s for at most one s′, P0 e_v = 0 otherwise;
  predicted dim = #arrows among {s, s′} other than the loop at s (0 if no s′); otherwise predicted 0."""
import sys, itertools
from cfibre import fibre
from extend import vertex_folds
def predicted(nv, arrows, P0):
    """generalized loop-collapse with spectators: for each looped s with P0 e_s = e_s whose row s is zero except at s
    and at collapsed vertices s′ with P0 e_{s′} = −e_s (column supported on s only); dim = Σ #arrows among {s, s′} other than the loop."""
    cols = [list(P0.column(v)) for v in range(nv)]; total = 0
    for s in range(nv):
        es = [1 if i == s else 0 for i in range(nv)]
        if cols[s] != es or (s, s) not in arrows: continue
        coll = [v for v in range(nv) if v != s and cols[v] == [-x for x in es]]
        other = [v for v in range(nv) if v != s and v not in coll and cols[v][s] != 0]
        if other or len(coll) != 1: continue
        sp = coll[0]; total += sum(1 for a in arrows if set(a) <= {s, sp} and a != (s, s))
    return total
def iso_classes(nv, loops_required=True):
    cand = [(i, j) for i in range(nv) for j in range(nv)]
    seen = set(); out = []
    for mask in range(1, 1 << len(cand)):
        E = [cand[k] for k in range(len(cand)) if mask >> k & 1]
        if loops_required and not any(a == b for a, b in E): continue
        canon = min(tuple(sorted((p[a], p[b]) for a, b in E)) for p in itertools.permutations(range(nv)))
        if canon not in seen: seen.add(canon); out.append(list(canon))
    return out
if __name__ == '__main__':
    nv, L = int(sys.argv[1]), int(sys.argv[2])
    qs = iso_classes(nv) if len(sys.argv) < 4 else [eval(x) for x in sys.argv[3].split(';')]
    allok = True
    for arrows in qs:
        A = set(arrows); rows = []
        for P0 in vertex_folds(nv):
            if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                       all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
            R, J, n = fibre(nv, arrows, P0, L)
            dim = J.dimension(); pr = predicted(nv, arrows, P0)
            rows.append((dim, pr, [list(P0.column(v)) for v in range(nv)]))
        ok = all(d == p for d, p, _ in rows); allok &= ok
        mism = [(d, p, c) for d, p, c in rows if d != p]
        print(f"{arrows}: {len(rows)} vertex folds; positive fibres {sorted(d for d, p, _ in rows if d > 0)}; prediction matches: {ok}{'' if ok else ' MISMATCH ' + str(mism)}", flush=True)
    print(f"collapse classification matches on all {len(qs)} quivers: {allok}")

"""Parallel, resumable census of Δ-compatible folds of kQ (Q acyclic, any parallel classes) over F_p.

Same mathematics as gcensus.py (kernel fibration; one kernel fold per grading-torus orbit; K-components enumerated one
column at a time with progressively filtered masked checks; full checks at the last column).
Parallel structure: a task is (fibre, partial fold X, column). A worker extends X by one column and returns the surviving
children, so one huge fibre is split across all workers. Tasks are scheduled depth-first so fibres finish early.
Checkpoint: completed fibres are saved (atomically, at most every 30 s); a rerun resumes from them.
Output: one JSON line with N, fibre histogram, orbit counts/sizes under G (column-sum-1 substitutions on every parallel
class) and under G × grading.  Usage:
  python3 pcensus.py QUIVER p [--jobs J] [--checkpoint FILE] [--out results.jsonl]"""
import sys, os, json, time, pickle, argparse, itertools, collections
import numpy as np
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gcensus.py')).read().split("def census(")[0])
QUIVERS.update({
    'K3':        (2, [(0, 1)] * 3),                                     # 1⇛2  (v2.18: 34 orbits)
    'K3_after':  (3, [(0, 1)] * 3 + [(1, 2)]),                          # 1⇛2→3
    'K3_srcout': (3, [(0, 1)] * 3 + [(0, 2)]),                          # 1⇛2, 1→3
    'two_pairs': (3, [(0, 1), (0, 1), (1, 2), (1, 2)]),                 # 1⇉2⇉3
    'two_pairs_sep': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3)]),     # 1⇉2→3⇉4
    'pair_after3': (5, [(0, 1), (0, 1), (1, 2), (2, 3), (3, 4)]),       # 1⇉2→3→4→5
})

C = {}   # worker context

def setup(name, p):
    nv, arrows = QUIVERS[name]; F = Field(p)
    basis, M, DL = build(nv, arrows); d = len(basis); I = np.eye(d, dtype=np.int64)
    Dmat = DL.reshape(d, d * d).T % p
    Kb = F.nullspace(Dmat); kdim = len(Kb)
    Cidx = []; cur = list(Kb)
    for i in range(d):
        if F.rank(np.array(cur + [I[i]])) > len(cur): cur.append(I[i]); Cidx.append(i)
    B = np.array(cur).T % p; Binv = F.inv(B)
    DC = Dmat[:, Cidx]; assert F.rank(DC) == len(Cidx)
    Kmat = np.array(Kb).T
    allK = np.array(list(itertools.product(range(p), repeat=kdim)), dtype=np.int64) @ np.array(Kb) % p
    return dict(name=name, p=p, F=F, basis=basis, arrows=arrows, M=M, DL=DL, d=d, I=I, Kb=Kb, kdim=kdim, Cidx=Cidx,
                B=B, Binv=Binv, DC=DC, Kmat=Kmat, allK=allK)

def init_worker(name, p):
    C.update(setup(name, p))

def kernel_folds_staged(ctx):
    """(Slower alternative, kept for cross-checks.) Algebra folds of K = ker Δ, enumerated column by column (K coordinates), with the same masked checks as the fibres."""
    p, F, kdim, Kmat, M, d = ctx['p'], ctx['F'], ctx['kdim'], ctx['Kmat'], ctx['M'], ctx['d']
    Kinv = ctx['Binv'][:kdim]                                            # old coords (on K) -> K coords
    MK = np.zeros((kdim, kdim, kdim), dtype=np.int64)
    for i in range(kdim):
        for j in range(kdim):
            prod = np.einsum('i,j,ijk->k', Kmat[:, i], Kmat[:, j], M) % p
            MK[i, j] = Kinv @ prod % p
            assert ((Kmat @ MK[i, j] - prod) % p == 0).all(), "K is not closed under multiplication"
    I = np.eye(kdim, dtype=np.int64)
    allv = np.array(list(itertools.product(range(p), repeat=kdim)), dtype=np.int64)
    def kcheck(Xb, known, new):
        st = {'Xb': Xb, 'idx': np.arange(len(Xb))}
        def Pap(v):
            Xb = st['Xb']; unk = v[..., known:].any(axis=-1)
            cc = v[:, :known] if v.ndim == 2 else np.broadcast_to(v[:known], (Xb.shape[0], known))
            return np.einsum('nij,nj->ni', Xb[:, :, :known], cc) % p, unk
        def keep(ok): st['Xb'] = st['Xb'][ok]; st['idx'] = st['idx'][ok]; return len(st['idx'])
        col = lambda i: st['Xb'][:, :, i] if i == new else st['Xb'][0, :, i]
        for i in range(known):
            v = col(i); r, unk = Pap(v)
            if not keep(~(((r - v) % p).any(axis=1) & ~unk)): return st['idx']
        for i in range(known):
            for j in range(known):
                for comp in (False, True):
                    xi, yj = col(i), col(j)
                    if comp: xi, yj = (I[:, i] - xi) % p, (I[:, j] - yj) % p
                    if xi.ndim == 1 and yj.ndim == 1: h = np.einsum('i,j,ijk->k', xi, yj, MK) % p
                    elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, MK) % p
                    elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, MK) % p
                    else: h = np.einsum('ni,nj,ijk->nk', xi, yj, MK) % p
                    r, unk = Pap(h)
                    if not keep(~((r if comp else (r - h) % p).any(axis=1) & ~unk)): return st['idx']
        return st['idx']
    layer = [np.zeros((kdim, kdim), dtype=np.int64)]
    for c in range(kdim):
        nxt = []
        for X in layer:
            Xb = np.repeat(X[None], len(allv), axis=0); Xb[:, :, c] = allv
            nxt += list(Xb[kcheck(Xb, c + 1, c)])
        layer = nxt
    return layer

def kernel_folds(ctx):
    """Algebra folds of K = ker Δ as complementary pairs of subalgebras (subspace enumeration)."""
    p, F, kdim, Kmat, M = ctx['p'], ctx['F'], ctx['kdim'], ctx['Kmat'], ctx['M']
    mul = lambda x, y: np.einsum('i,j,ijk->k', x, y, M) % p
    def subspaces(n, k):
        for piv in itertools.combinations(range(n), k):
            free = [(i, j) for i in range(k) for j in range(n) if j not in piv and j > piv[i]]
            for vals in itertools.product(range(p), repeat=len(free)):
                S = np.zeros((k, n), dtype=np.int64)
                for i, c in enumerate(piv): S[i, c] = 1
                for (i, j), v in zip(free, vals): S[i, j] = v
                yield S
    def is_subalg(S):
        vecs = [Kmat @ r % p for r in S]
        return all(F.rank(np.array(vecs + [mul(x, y)])) == len(vecs) for x in vecs for y in vecs) if vecs else True
    subalgs = [S for k in range(kdim + 1) for S in subspaces(kdim, k) if is_subalg(S)]
    out = []
    for H in subalgs:
        for Q in subalgs:
            if H.shape[0] + Q.shape[0] != kdim: continue
            Bm = np.vstack([H, Q]).reshape(kdim, kdim)
            if F.rank(Bm) != kdim: continue
            Bt = Bm.T % p; E = np.diag([1] * H.shape[0] + [0] * Q.shape[0])
            out.append(Bt @ E @ F.inv(Bt) % p)
    return out

def lift(X, col):
    p, d, I, DL, Binv, DC, Cidx, kdim, F = (C[k] for k in ('p', 'd', 'I', 'DL', 'Binv', 'DC', 'Cidx', 'kdim', 'F'))
    T = DL[Cidx[col - kdim]]; acc = np.zeros((d, d), dtype=np.int64)
    for a, b in zip(*np.nonzero(T)):
        ca, cb = Binv @ I[a] % p, Binv @ I[b] % p
        assert not ca[col:].any() and not cb[col:].any(), "Δ-lift needs an unknown column"
        acc = (acc + T[a, b] * np.outer(X[:, :col] @ ca[:col] % p, X[:, :col] @ cb[:col] % p)) % p
    sol = F.solve(DC, acc.ravel())
    if sol is None: return None
    return sum((sol[k] * I[Cidx[k]] for k in range(len(Cidx))), np.zeros(d, dtype=np.int64)) % p

def check(Xb, known, new):
    p, M, Binv, B = C['p'], C['M'], C['Binv'], C['B']
    st = {'Xb': Xb, 'idx': np.arange(Xb.shape[0])}
    def Papply(v):
        Xb = st['Xb']; c = np.einsum('ij,...j->...i', Binv, v) % p; unk = c[..., known:].any(axis=-1)
        cc = c[:, :known] if c.ndim == 2 else np.broadcast_to(c[:known], (Xb.shape[0], known))
        return np.einsum('nij,nj->ni', Xb[:, :, :known], cc) % p, unk
    def keep(ok):
        st['Xb'] = st['Xb'][ok]; st['idx'] = st['idx'][ok]; return len(st['idx'])
    col = lambda i: st['Xb'][:, :, i] if i == new else st['Xb'][0, :, i]
    for i in range(known):
        v = col(i); r, unk = Papply(v)
        if not keep(~(((r - v) % p).any(axis=1) & ~unk)): return st['idx']
    for i in range(known):
        for j in range(known):
            for comp in (False, True):
                xi, yj = col(i), col(j)
                if comp: xi, yj = (B[:, i] - xi) % p, (B[:, j] - yj) % p
                if xi.ndim == 1 and yj.ndim == 1: h = np.einsum('i,j,ijk->k', xi, yj, M) % p
                elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, M) % p
                elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, M) % p
                else: h = np.einsum('ni,nj,ijk->nk', xi, yj, M) % p
                r, unk = Papply(h)
                if not keep(~((r if comp else (r - h) % p).any(axis=1) & ~unk)): return st['idx']
    return st['idx']

def extend(task):
    """task = (rep, X, col) -> (rep, col, children): partial folds with columns 0..col filled that pass the masked checks."""
    rep, X, col = task
    X = X.astype(np.int64); L = lift(X, col)
    if L is None: return rep, col, []
    Xb = np.repeat(X[None], len(C['allK']), axis=0); Xb[:, :, col] = (L + C['allK']) % C['p']
    kids = Xb[check(Xb, col + 1, col)]
    return rep, col, [k.astype(np.int16) for k in kids]

def orbit_label(folds, mats, p, chunk=4096):
    """G-orbits: dict fold-key -> orbit id, list of representatives, list of orbit sizes."""
    g, gi = mats; lab = {}; reps = []; sizes = []; key = lambda A: (A % p).astype(np.int8).tobytes()
    for P in folds:
        if key(P) in lab: continue
        o = set()
        for s_ in range(0, len(g), chunk):
            Q = np.einsum('gij,jk,gkl->gil', g[s_:s_ + chunk], P, gi[s_:s_ + chunk]) % p
            o.update(q.astype(np.int8).tobytes() for q in Q)
        for k in o: lab[k] = len(reps)
        reps.append(P); sizes.append(len(o))
    return lab, reps, sizes

def merge(lab, reps, sizes, extra, p):
    """orbits of ⟨G, extra⟩ (extra normalizes G): union of G-orbits under the extra elements. Returns count, size histogram,
    and one representative fold per merged orbit."""
    par = list(range(len(reps)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    key = lambda A: (A % p).astype(np.int8).tobytes()
    for i, P in enumerate(reps):
        for t, ti in extra:
            j = lab[key(t @ P @ ti)]; a, b = find(i), find(j)
            if a != b: par[a] = b
    comp = collections.defaultdict(int); rep = {}
    for i in range(len(reps)):
        r = find(i); comp[r] += sizes[i]; rep.setdefault(r, reps[i])
    return len(comp), dict(sorted(collections.Counter(comp.values()).items())), list(rep.values())

def orbit_stats(folds, mats, p, chunk=4096):
    """orbit count and size histogram of the group given as stacked (g, g⁻¹) arrays; vectorized over the group."""
    g, gi = mats; seen = set(); sizes = collections.Counter(); key = lambda A: (A % p).astype(np.int8).tobytes()
    for P in folds:
        if key(P) in seen: continue
        o = set()
        for s in range(0, len(g), chunk):
            Q = np.einsum('gij,jk,gkl->gil', g[s:s + chunk], P, gi[s:s + chunk]) % p
            o.update(q.astype(np.int8).tobytes() for q in Q)
        seen |= o; sizes[len(o)] += 1
    return sum(sizes.values()), dict(sorted(sizes.items()))

def lie_G(ctx):
    """Lie(G) on A: for each parallel class (a_1..a_m), derivations a_j ↦ a_i − a_m (i < m), extended to paths; plus E (grading)."""
    basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']; idxb = {b: i for i, b in enumerate(basis)}
    classes = collections.defaultdict(list)
    for ai, (s_, t_) in enumerate(arrows): classes[(s_, t_)].append(ai)
    gens = []
    for c in (c for c in classes.values() if len(c) > 1):
        for i in range(len(c) - 1):
            for j in range(len(c)):
                D = np.zeros((d, d), dtype=np.int64)          # derivation with D(a_j) = a_i − a_m
                for col, b in enumerate(basis):
                    if b[0] == 'v': continue
                    for pos, x in enumerate(b[1]):
                        if x != c[j]: continue
                        for y, coef in ((c[i], 1), (c[-1], -1)):
                            D[idxb[('p', b[1][:pos] + (y,) + b[1][pos + 1:])], col] += coef
                gens.append(D)
    E = np.diag([0 if b[0] == 'v' else len(b[1]) for b in basis]).astype(np.int64)
    return gens, E

def transversal(reps, ctx):
    """For one representative per (G × grading)-orbit: r_G = rank{[D,P]}, r_H = rank with E added (mod p, and over Q
    when the representative lifts to an exact Q-fold with entries p−1 ↦ −1). Counts orbits with r_H > r_G."""
    import sympy as sp
    p, F = ctx['p'], ctx['F']; gens, E = lie_G(ctx)
    cnt = collections.Counter()
    M, DL, d = ctx['M'], ctx['DL'], ctx['d']
    def q_fold(P):
        Pq = np.where(P == p - 1, -1, P).astype(np.int64)
        if (Pq @ Pq != Pq).any(): return None
        I = np.eye(d, dtype=np.int64)
        for i in range(d):
            for j in range(d):
                h = np.einsum('i,j,ijk->k', Pq[:, i], Pq[:, j], M); q = np.einsum('i,j,ijk->k', (I - Pq)[:, i], (I - Pq)[:, j], M)
                if (Pq @ h != h).any() or (Pq @ q).any(): return None
        for l in range(d):
            lhs = np.einsum('m,mab->ab', Pq[:, l], DL)
            if (lhs != Pq @ DL[l] @ Pq.T).any(): return None
        return Pq
    for P in reps:
        rows = [((D @ P - P @ D) % p).ravel() for D in gens]
        rG = F.rank(np.array(rows)) if rows else 0; rH = F.rank(np.array(rows + [((E @ P - P @ E) % p).ravel()]))
        Pq = q_fold(P); tag = 'modp'
        if Pq is not None:
            rowsq = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
            rGq = sp.Matrix(rowsq).rank() if rowsq else 0; rHq = sp.Matrix(rowsq + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
            cnt['lifted'] += 1; cnt['lifted_transversal'] += int(rHq > rGq)
        cnt['transversal_modp'] += int(rH > rG)
    return dict(cnt)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('p', type=int)
    ap.add_argument('--jobs', type=int, default=os.cpu_count()); ap.add_argument('--checkpoint'); ap.add_argument('--out')
    ap.add_argument('--no-orbits', action='store_true'); ap.add_argument('--save-folds')
    a = ap.parse_args(); name, p = a.quiver, a.p
    t0 = time.time(); ctx = setup(name, p); C.update(ctx)
    d, kdim, B, Binv, Kmat, basis = ctx['d'], ctx['kdim'], ctx['B'], ctx['Binv'], ctx['Kmat'], ctx['basis']
    Kf = kernel_folds(ctx)
    lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
    phis = {l: (np.diag([pow(l, k, p) for k in lens]), np.diag([pow(pow(l, p - 2, p), k, p) for k in lens])) for l in range(1, p)}
    kkey = lambda A: (A % p).tobytes(); done_k = set(); reps = []      # (PK, {kkey: l})
    for PK in Kf:
        if kkey(PK) in done_k: continue
        torb = {}
        for l, (ph, phi_) in phis.items():
            PKl = (Binv @ ph @ (Kmat @ PK % p) % p)[:kdim] % p
            PKl = PKl @ (Binv @ phi_ @ Kmat % p)[:kdim] % p
            torb.setdefault(kkey(PKl), l)
        done_k |= set(torb); reps.append((PK, list(torb.values())))
    print(f"{name} F_{p}: d={d}, dim K={kdim}, {len(Kf)} kernel folds in {len(reps)} torus orbits; jobs={a.jobs}", flush=True)
    fibres = {}
    if a.checkpoint and os.path.exists(a.checkpoint):
        ck = pickle.load(open(a.checkpoint, 'rb'))
        assert (ck['name'], ck['p'], ck['nreps']) == (name, p, len(reps)), "checkpoint is for another run"
        fibres = ck['fibres']; print(f"  resumed: {len(fibres)}/{len(reps)} fibres done", flush=True)
    last_save = [time.time()]
    def save(force=False):
        if not a.checkpoint or (not force and time.time() - last_save[0] < 30): return
        tmp = a.checkpoint + '.tmp'
        pickle.dump(dict(name=name, p=p, nreps=len(reps), fibres=fibres), open(tmp, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, a.checkpoint); last_save[0] = time.time()
    stack = []; outstanding = collections.Counter(); partial = collections.defaultdict(list)
    todo = [r for r in range(len(reps)) if r not in fibres]
    for r in reversed(todo):
        X = np.zeros((d, d), dtype=np.int64); X[:, :kdim] = Kmat @ reps[r][0] % p
        if kdim == d: fibres[r] = [X.astype(np.int16)]; continue
        stack.append((r, X.astype(np.int16), kdim)); outstanding[r] += 1
    with ProcessPoolExecutor(a.jobs, initializer=init_worker, initargs=(name, p)) as ex:
        inflight = set(); lastlog = time.time(); ntasks = 0
        while stack or inflight:
            while stack and len(inflight) < 4 * a.jobs:
                inflight.add(ex.submit(extend, stack.pop()))
            fin, inflight = wait(inflight, return_when=FIRST_COMPLETED)
            for fu in fin:
                r, col, kids = fu.result(); outstanding[r] -= 1; ntasks += 1
                if col + 1 == d: partial[r] += kids
                else:
                    for k in kids: stack.append((r, k, col + 1)); outstanding[r] += 1
                if outstanding[r] == 0:
                    fibres[r] = partial.pop(r, []); del outstanding[r]; save()
            if time.time() - lastlog > 60:
                print(f"  {len(fibres)}/{len(reps)} fibres done, {len(outstanding)} open, {len(stack)} queued tasks, "
                      f"{ntasks} tasks run, t={time.time()-t0:.0f}s", flush=True); lastlog = time.time()
    save(force=True)
    # materialize all folds (old coordinates) via the grading on each torus orbit
    folds = []; fib = collections.Counter()
    for r, (PK, ls) in enumerate(reps):
        base = [X.astype(np.int64) @ Binv % p for X in fibres[r]]
        for l in ls:
            ph, phi_ = phis[l]; fib[len(base)] += 1
            folds += [ph @ P_ @ phi_ % p for P_ in base]
    assert sum(fib.values()) == len(Kf)
    assert len({P_.astype(np.int8).tobytes() for P_ in folds}) == len(folds)
    res = dict(name=name, p=p, d=d, kdim=kdim, nK=len(Kf), reps=len(reps), N=len(folds), fib=dict(sorted(fib.items())))
    if not a.no_orbits:
        # G: column-sum-1 substitutions on every parallel class (as in gcensus.py), stacked
        F = ctx['F']; arrows = ctx['arrows']
        classes = collections.defaultdict(list)
        for ai, (s_, t_) in enumerate(arrows): classes[(s_, t_)].append(ai)
        classes = [c for c in classes.values() if len(c) > 1]
        def mats(m):      # column sums 1: the last row is determined by the first m − 1
            out = []
            for vals in itertools.product(range(p), repeat=(m - 1) * m):
                top = np.array(vals, dtype=np.int64).reshape(m - 1, m)
                g = np.vstack([top, (1 - top.sum(axis=0)) % p])
                if F.rank(g) == m: out.append(g)
            return out
        idxb = {b: i for i, b in enumerate(basis)}
        def phi(subs):
            Phi = np.zeros((d, d), dtype=np.int64)
            for i, b in enumerate(basis):
                if b[0] == 'v': Phi[i, i] = 1; continue
                for choice in itertools.product(*[list(subs.get(x, {x: 1}).items()) for x in b[1]]):
                    j = idxb[('p', tuple(x for x, _ in choice))]
                    Phi[j, i] = (Phi[j, i] + int(np.prod([c for _, c in choice]))) % p
            return Phi
        G = []
        for gs in itertools.product(*[mats(len(c)) for c in classes]):
            subs = {}
            for c, g in zip(classes, gs):
                for jc, x in enumerate(c): subs[x] = {c[ir]: int(g[ir, jc]) for ir in range(len(c)) if g[ir, jc]}
            Phi = phi(subs); G.append((Phi, F.inv(Phi)))
        Gg = (np.array([g for g, _ in G]), np.array([gi for _, gi in G]))
        lab, greps, gsizes = orbit_label(folds, Gg, p)
        res['G'] = len(G); res['orbits_G'] = len(greps); res['sizes_G'] = dict(sorted(collections.Counter(gsizes).items()))
        res['orbits_GxGr'], res['sizes_GxGr'], hreps = merge(lab, greps, gsizes, list(phis.values()), p)
        res['transversal'] = transversal(hreps, ctx)
        # class-weight torus: φ_w scales each parallel class (single arrows included) by w_c; keep the w preserving Folds
        cls = collections.defaultdict(list)
        for ai, st in enumerate(arrows): cls[st].append(ai)
        cl = list(cls.values()); FS = set(lab); Tpres = []
        for w in itertools.product(range(1, p), repeat=len(cl)):
            wa = {a: w[k] for k, c in enumerate(cl) for a in c}
            diag = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
            ph, phi_ = np.diag(diag), np.diag([pow(x, p - 2, p) for x in diag])
            if all(((ph @ P_ @ phi_) % p).astype(np.int8).tobytes() in FS for P_ in folds): Tpres.append((ph, phi_))
        res['classes'] = len(cl); res['torus_preserving'] = len(Tpres); res['torus_all'] = (p - 1) ** len(cl)
        res['orbits_GxT'], res['sizes_GxT'], _ = merge(lab, greps, gsizes, Tpres, p)
    res['t'] = round(time.time() - t0)
    line = json.dumps(res); print(line, flush=True)
    if a.out: open(a.out, 'a').write(line + "\n")
    if a.save_folds: pickle.dump(dict(name=name, p=p, folds=[P_.astype(np.int8) for P_ in folds]), open(a.save_folds, 'wb'))

if __name__ == '__main__':
    main()

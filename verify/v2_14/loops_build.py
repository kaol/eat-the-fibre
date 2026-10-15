import numpy as np
def tbuild(nv, arrows, LMAX):
    allp = [(v,) for v in range(nv)]; frontier = [('a', k) for k in range(len(arrows))]; L = 1
    while frontier and L <= LMAX:
        allp += frontier; nf = []
        for p in frontier:
            last = arrows[p[-1]][1]
            for k, (s, t) in enumerate(arrows):
                if s == last: nf.append(p + (k,))
        frontier = nf; L += 1
    idx = {p: i for i, p in enumerate(allp)}; d = len(allp)
    start = lambda p: p[0] if p[0] != 'a' else arrows[p[1]][0]
    end = lambda p: p[0] if p[0] != 'a' else arrows[p[-1]][1]
    mult = np.zeros((d, d, d)); Dl = np.zeros((d, d, d))
    for p in allp:
        for q in allp:
            if end(p) != start(q): continue
            r = q if p[0] != 'a' else (p if q[0] != 'a' else p + q[1:])
            if r in idx: mult[idx[p], idx[q], idx[r]] = 1
    for p in allp:
        if p[0] != 'a': continue
        ks = p[1:]
        for i, k in enumerate(ks):
            pre = ('a',) + ks[:i] if i > 0 else (arrows[k][0],)
            suf = ('a',) + ks[i + 1:] if i + 1 < len(ks) else (arrows[k][1],)
            Dl[idx[p], idx[pre], idx[suf]] += 1
    return d, mult, Dl, allp

Q = {'loop': (1, [(0, 0)]), 'loop+arrow': (2, [(0, 0), (0, 1)]), 'arrow+loop at target': (2, [(0, 1), (1, 1)]),
     'two loops+arrow': (2, [(0, 0), (0, 1), (1, 1)]), '2-cycle': (2, [(0, 1), (1, 0)]), '3-cycle': (3, [(0, 1), (1, 2), (2, 0)]),
     '2-cycle+loop': (2, [(0, 1), (1, 0), (0, 0)]), '2-cycle+loop+tail': (3, [(0, 1), (1, 0), (0, 0), (1, 2), (2, 0)]), 'triangle+loop': (3, [(0, 1), (1, 2), (2, 0), (0, 0), (1, 0)]), 'complete digraph K2 with loops': (2, [(0, 0), (0, 1), (1, 0), (1, 1)])}

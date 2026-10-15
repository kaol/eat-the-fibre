"""Polynomial-time count of arrow-compatible splittings via the leaf-type rule (v2.6 Thm 3)."""
import sys, functools
from math import comb, factorial
from collections import defaultdict
@functools.lru_cache(None)
def f(x, y):   # bipartite rooted forests on x + y vertices = spanning trees of K_{1,x,y}
    if x == 0 and y == 0: return 1
    n = 1 + x + y
    from fractions import Fraction as Fr
    v = Fr(n) * Fr(1 + y) ** (x - 1) * Fr(1 + x) ** (y - 1)
    assert v.denominator == 1; return int(v)
@functools.lru_cache(None)
def Fleaf(bh, bp, X, Y):   # forests on B where specified X (H) and Y (P) vertices have no B-children
    return f(bh - X, bp - Y) * (bp - Y + 1) ** X * (bh - X + 1) ** Y
@functools.lru_cache(None)
def G(bh, bp, zH, zP):     # forests on B with exactly zH H-vertices and zP P-vertices childless (within B)
    s = 0
    for i in range(bh - zH + 1):
        for j in range(bp - zP + 1):
            s += (-1) ** (i + j) * comb(bh - zH, i) * comb(bp - zP, j) * Fleaf(bh, bp, zH + i, zP + j)
    return comb(bh, zH) * comb(bp, zP) * s
@functools.lru_cache(None)
def AttEP(items, bh, zH):  # items -> B_H ∪ root, zH specified B_H bins nonempty
    return sum((-1) ** i * comb(zH, i) * (bh + 1 - i) ** items for i in range(zH + 1))
@functools.lru_cache(None)
def AttHL(items, bp, zP):  # items -> B_P ∪ root, zP specified B_P bins get >= 2 items
    s = 0
    for i in range(zP + 1):
        for l in range(zP - i + 1):
            if l > items: continue
            ff = factorial(items) // factorial(items - l)
            s += (-1) ** (i + l) * comb(zP, i) * comb(zP - i, l) * ff * (bp + 1 - i - l) ** (items - l)
    return s
@functools.lru_cache(None)
def Nforest(e, h, m, bh, bp):
    tot = 0
    for zH in range(bh + 1):
        for zP in range(bp + 1):
            g = G(bh, bp, zH, zP)
            if g: tot += g * AttEP(e + m, bh, zH) * AttHL(h - m, bp, zP)
    return tot
@functools.lru_cache(None)
def Ntotal(e, h, m, b):
    return sum(comb(b, bh) * Nforest(e, h, m, bh, b - bh) for bh in range(b + 1))
def a(n):
    # DP over positions. symbols: 'E', 'B', 'H' (unmarked HL), 'K' (marked HL block element), 'P' (P1 run element)
    # state: (prev, curP, curK, unmatched_runs tuple, unmatched_blocks tuple, e, h, m, b)
    # curP / curK = length of the current P1-run / marked block (0 if not inside one)
    from collections import Counter
    start = ('S', 0, 0, (), (), 0, 0, 0, 0)
    dp = {start: 1}
    def close(lst_runs, lst_blocks, kind, L):
        """close a run (kind 'P') or block (kind 'K') of length L: match with an unmatched opposite item, or leave unmatched."""
        out = []
        runs, blocks = list(lst_runs), list(lst_blocks)
        if kind == 'P':
            c = blocks.count(L)
            if c:
                bl = list(blocks); bl.remove(L); out.append((tuple(runs), tuple(sorted(bl)), c))
            out.append((tuple(sorted(runs + [L])), tuple(blocks), 1))
        else:
            c = runs.count(L)
            if c:
                rl = list(runs); rl.remove(L); out.append((tuple(sorted(rl)), tuple(blocks), c))
            out.append((tuple(runs), tuple(sorted(blocks + [L])), 1))
        return out
    for pos in range(n):
        nd = defaultdict(int)
        for (prev, cP, cK, R, Bl, e, h, m, b), w in dp.items():
            for sym in 'EBHKP':
                # adjacency: consecutive non-E pairs allowed only (HL-type,HL-type) or (P,P) in the same run
                hl = lambda s: s in 'HK'
                if prev not in ('S', 'E') and sym != 'E':
                    if not ((hl(prev) and hl(sym)) or (prev == 'P' and sym == 'P')): continue
                # closing current run/block when it does not continue
                opts = [(R, Bl, 1)]
                ncP, ncK = 0, 0
                if prev == 'P':
                    if sym == 'P': ncP = cP + 1; opts = [(R, Bl, 1)]
                    else: opts = close(R, Bl, 'P', cP)
                if prev == 'K':
                    if sym == 'K':
                        # either continue the block or close it and start a new marked block
                        for (R2, B2, c2) in [(R, Bl, 1)]:
                            pass
                    pass
                if sym == 'P' and prev != 'P': ncP = 1
                for (R2, B2, c2) in opts:
                    if prev == 'K':
                        if sym == 'K':
                            # continue block
                            st = (sym, ncP, cK + 1, R2, B2, e, h + 1, m, b); nd[st] += w * c2
                            # or close current block and start a new one
                            for (R3, B3, c3) in close(R2, B2, 'K', cK):
                                st = (sym, ncP, 1, R3, B3, e, h + 1, m, b); nd[st] += w * c2 * c3
                            continue
                        else:
                            for (R3, B3, c3) in close(R2, B2, 'K', cK):
                                st = (sym, ncP, 1 if sym == 'K' else 0, R3, B3,
                                      e + (sym == 'E'), h + (sym in 'HK'), m + (sym == 'P'), b + (sym == 'B'))
                                nd[st] += w * c2 * c3
                            continue
                    st = (sym, ncP, 1 if sym == 'K' else 0, R2, B2,
                          e + (sym == 'E'), h + (sym in 'HK'), m + (sym == 'P'), b + (sym == 'B'))
                    nd[st] += w * c2
        dp = nd
    total = 0
    for (prev, cP, cK, R, Bl, e, h, m, b), w in dp.items():
        opts = [(R, Bl, 1)]
        if prev == 'P': opts = [(r, bb, c) for (r, bb, c) in close(R, Bl, 'P', cP)]
        if prev == 'K': opts = [(r, bb, c) for (r, bb, c) in close(R, Bl, 'K', cK)]
        for (R2, B2, c2) in opts:
            if R2 == () and B2 == (): total += w * c2 * Ntotal(e, h, m, b)
    return total
if __name__ == '__main__':
    known = [2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445, 25665390]
    vals = [a(n) for n in range(1, int(sys.argv[1]) + 1)]
    print(vals)
    print("matches brute force n ≤ 10:", vals[:10] == known[:len(vals[:10])])

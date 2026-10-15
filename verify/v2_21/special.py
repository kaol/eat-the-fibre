"""Geometry of the special fibre (over the kernel fold P_K = id on e1,e2,e3, δ ↦ 0) of 1⇉2→3 over F_p:
G = Aff(1) fixes P_K, so it acts on the fibre; orbit decomposition 4 free + 1 of size p accounts for 4p² − 3p."""
import sys, itertools, collections
exec(open("q123s2.py").read().split("exec(\"t0 = time.time(); done\"")[0])
PK = np.diag([1, 1, 1, 0]).astype(np.int64)
F = fibre(PK); S = {f.tobytes() for f in F}; X = [f.ravel() for f in F]
print(f"p={p}: special fibre has {len(F)} points (4p²−3p = {4*p*p-3*p})")
pairs = [(L8(B), L8(inv[B.tobytes()])) for B in Gb]
seen = set(); osz = collections.Counter()
for f in F:
    if f.tobytes() in seen: continue
    o = {((g @ f @ gi) % p).tobytes() for g, gi in pairs}; seen |= o; osz[len(o)] += 1
print(f"G-orbits on the special fibre: sizes {dict(sorted(osz.items()))}  (4 free + 1 of size p predicts {{{p}: 1, {p*(p-1)}: 4}})")

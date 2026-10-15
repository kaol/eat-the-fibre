# Orbit-size distribution of G(F_p) on Δ-compatible folds of k(1⇉2→3); reuses q123v census.
import sys, pickle, collections
src = open('q123v.py').read()
exec(src.split('exec("G = []"')[0])
exec("G = []" + open('q123.py').read().split("G = []")[1].split("seen = set()")[0])
seen = set(); sizes = collections.Counter(); reps = []
for P in folds:
    if P.tobytes() in seen: continue
    orb = {((g @ P @ gi) % p).tobytes() for g, gi in zip(G, Ginv)}
    seen |= orb; sizes[len(orb)] += 1; reps.append(P)
print(f"p={p} |G|={len(G)} orbits={sum(sizes.values())} folds={len(folds)} size distribution={dict(sorted(sizes.items()))}")
pickle.dump(dict(p=p, folds=folds, reps=reps), open(f'folds_p{p}.pkl', 'wb'))

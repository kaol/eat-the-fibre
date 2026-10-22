"""v2.45: local dimension of behaviours with fixed image, at lifted representatives of given signature groups.

For a fold P₀ (integer lift of a census representative) with image A = im P₀:
  W = {P₀ Y (1 − P₀)} = Hom(ker P₀, A)            (tilts of the kernel; P₀ + N is idempotent with image A)
  T_V = W ∩ ker J                                   (tangent of V_A, the folds with image A; J = Jacobian of Folds)
  Der_A = {D ∈ Der(kQ) : (1 − P₀) D P₀ = 0}         (derivations preserving A)
  m = dim T_V − dim([Der_A, P₀] ∩ T_V)
and v2.34's b = dim T_P Folds − dim([Der, P₀] ∩ T). Linear algebra modulo a 7-digit prime (ranks agree with ℚ for all
but finitely many primes); the census-prime (b, e) of the orbit is reported alongside.
Usage: python3 analysis/v2_45_fixed_image.py REPS_p.pkl SIG_SUBSTRING [MAX_PER_CLASS]"""
import os, sys, io, contextlib, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rf, tg = sys.argv[1], sys.argv[2]; maxpc = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
PR = 9999991
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
r, sig, p, d = g['r'], g['sig_of_loc'], g['p'], g['d']
M = np.array(g['M'], np.int64); DL = np.array(g['DL'], np.int64)
Cp = fr['setup'](r['name'], p); CP = fr['setup'](r['name'], PR)
I = np.eye(d, dtype=np.int64)
rank, nullspace = fr['rank'], fr['nullspace']
def inter_dim(A, B):
    """dim(rowspace A ∩ rowspace B)"""
    if not len(A) or not len(B): return 0
    return rank(A, PR) + rank(B, PR) - rank(np.concatenate([A, B]), PR)
out = collections.defaultdict(list); per = collections.Counter()
done = set()
for i, l in enumerate(r['loc']):
    if l in done or tg not in str(sig[l]): continue
    Pz = g['lifted'](i)
    if Pz is None: continue
    done.add(l)
    a = fr['at'](r['greps'][i], Cp, p); be = (a['b'], a['e'])
    if per[be] >= maxpc: continue
    per[be] += 1
    P0 = Pz.astype(np.int64) % PR; Q0 = (I - P0) % PR
    J = fr['jacobian'](P0, dict(d=d, M=M, DL=DL), PR)
    # W basis: P0 E_ij Q0 for elementary E_ij, then reduce to a basis
    Wv = np.array([(P0[:, [a_]] @ Q0[[b_], :] % PR).ravel() for a_ in range(d) for b_ in range(d)])
    Wb = Wv[[k for k in range(len(Wv)) if Wv[k].any()]]
    R, piv = fr['rref'](Wb, PR); Wb = R[:len(piv)]                      # basis of W
    # T_V = W ∩ ker J: coefficients c with J (c·Wb)ᵀ = 0
    Cc = nullspace((J @ Wb.T) % PR, PR)
    TV = (Cc @ Wb) % PR if len(Cc) else np.zeros((0, d * d), np.int64)
    # Der_A: D with Q0 D P0 = 0
    Der = CP['Der']
    cond = np.array([((Q0 @ D @ P0) % PR).ravel() for D in Der]).T
    cA = nullspace(cond % PR, PR)
    DerA = [np.einsum('k,kij->ij', c, Der) % PR for c in cA]
    brA = np.array([((D @ P0 - P0 @ D) % PR).ravel() for D in DerA]) if DerA else np.zeros((0, d * d), np.int64)
    m = len(TV) - inter_dim(TV, brA)
    out[be].append(m)
    print(f"orbit {l} rep {i} (b, e) = {be}: dim W = {len(Wb)}, dim T_V = {len(TV)}, dim Der_A = {len(DerA)}, m = {m}", flush=True)
print("summary: (b, e) -> m values:", {k: sorted(collections.Counter(v).items()) for k, v in sorted(out.items())})

"""v2.61: the λ = 1 split locus W_V of Q₁ over each stratum S_V of Q₀, as a union of components with their fibre dimensions.

Reuses the Singular scripts written by v2_61_mu_linear.py (OUT_PREFIX_BLOCK.sing): after the saturation loop L (the
fold ideal, generated in μ-degree ≤ 1) is decomposed (minAssGTZ over ℚ, ring with vs and μ); for each component W the
dimension, the dimension of its projection to the stratum (elimination of μ), and the fibre dimension (difference) and
whether W contains the identity lift locus (μ = 0 on all of W).
Usage: python3 analysis/v2_61_mu_bundle.py OUT_PREFIX"""
import sys, os, glob, subprocess, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import singular_nohilb; singular_nohilb.install()          # minAssGTZ without the overflowing Hilbert series
pref = sys.argv[1]
for fn in sorted(glob.glob(pref + '_*.sing')):
    if fn.endswith('_bundle.sing'): continue
    src = open(fn).read()
    head = src.split('print("ROUNDS "')[0]
    mus = re.search(r"ring r = 0,\((.*)\),dp", src).group(1).split(',')
    mu = [v for v in mus if re.fullmatch(r"m\d", v)]
    tail = ['if (size(SL) == 1 && deg(SL[1]) == 0) { print("EMPTY"); quit; }',
            'LIB "primdec.lib"; list LL = minAssGTZ(SL); int jj;',
            'for (jj = 1; jj <= size(LL); jj++) { ideal W = std(LL[jj]); ideal E = std(eliminate(W, ' + '*'.join(mu) + '));',
            '  ideal Z0 = std(W + ideal(' + ','.join(mu) + ')); print("COMP " + string(jj) + " DIM " + string(dim(W) - 1) + " BASE " + string(dim(E) - 1 - ' + str(len(mu)) + ') + " MUZERO " + string(dim(Z0) == dim(W))); kill W; kill E; kill Z0; }',
            'quit;']
    out_fn = fn.replace('.sing', '_bundle.sing'); open(out_fn, 'w').write(head + "\n".join(tail))
    try: res = subprocess.run(["Singular", "-q", "--no-rc", out_fn], capture_output=True, text=True, timeout=3000)
    except subprocess.TimeoutExpired: print(fn, "timeout"); continue
    if "EMPTY" in res.stdout:
        blk = re.search(r"_(\d{9})\.sing$", fn).group(1)
        print(f"{blk[:3]},{blk[3:6]},{blk[6:]}: no λ = 1 split lift (the conditions are inconsistent: 1 ∈ L)", flush=True); continue
    rows = re.findall(r"COMP (\d+) DIM (-?\d+) BASE (-?\d+) MUZERO (\d)", res.stdout)
    blk = re.search(r"_(\d{9})\.sing$", fn).group(1)
    fib = [int(d) - int(b) for _, d, b, _ in rows]
    print(f"{blk[:3]},{blk[3:6]},{blk[6:]}: {len(rows)} components of the λ = 1 split locus; fibre dimensions {fib}; "
          f"identity-lift components {sum(1 for r in rows if r[3] == '1')}" + (" ERROR " + (res.stdout + res.stderr)[-200:] if "? " in res.stdout + res.stderr else ""), flush=True)

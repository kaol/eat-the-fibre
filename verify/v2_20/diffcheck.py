import sys, numpy as np
sys.argv = ['x', '3']
src = open('q123v.py').read().replace('exec("G = []" + open(\'q123.py\').read().split("G = []")[1])', '')
exec(src)
vec = {f.tobytes(): f for f in folds}
src2 = open('q123.py').read().split("G = []")[0].replace("LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else None", "LIMIT = None")
ns = {}; exec(src2, ns)
orig = {f.tobytes(): f for f in ns['folds']}
extra = [vec[k] for k in vec if k not in orig]; missing = [orig[k] for k in orig if k not in vec]
print("vectorized", len(vec), "original", len(orig), "extra in vectorized:", len(extra), "missing:", len(missing))
for f in extra:
    print("extra passes full_check (incl. Δ):", ns['full_check'](f))
    print(f)

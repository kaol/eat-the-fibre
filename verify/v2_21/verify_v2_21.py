"""verify_v2_21.py — machine checks W1–W9 for discordianism-open-questions-explored-v2.21.md (quiver 1⇉2→3)
Install:  pip install -r requirements.txt.
Default (a few minutes): W1–W7.  SLOW=1 adds the F_7 census and fibre histogram (W8) and the special-fibre orbits at p = 7.
SLOW=2 also adds the F_11 census (W9)."""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = int(os.environ.get('SLOW', '0') or 0)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(*args):
    r = subprocess.run([sys.executable, '-u', *args], capture_output=True, text=True)
    if r.returncode: print(f"  !! {' '.join(args)} failed:", *r.stderr.strip().splitlines()[-3:], sep="\n     ")
    return r.stdout
def parse(out):
    m = re.search(r'F_(\d+): (\d+) folds .*?orbits = (\d+), orbit sizes (\{[^}]*\}), non-length-filtered folds = (\d+)', out)
    h = re.search(r'fibre sizes over the \d+ kernel folds: (\{[^}]*\})', out)
    return dict(N=int(m.group(2)), orbits=int(m.group(3)), sizes=eval(m.group(4)), nonfilt=int(m.group(5)),
                fib=eval(h.group(1)) if h else None, line=m.group(0)) if m else None
def pred(q):
    return dict(N=21 * q * q + 6 * q + 12, orbits=71, sizes={1: 23, q - 1: 11, q: 16, q * (q - 1): 21}, nonfilt=q - 1,
                fib={0: 24 * q * q + 19 * q + 18, 1: 11 * q + 12, q: 16 * q - 2, q * q: 1, 4 * q * q - 3 * q: 1})
def matches(s, q):
    return s is not None and all(s[k] == pred(q)[k] for k in ('N', 'orbits', 'sizes', 'nonfilt', 'fib'))

s3 = parse(run('q123f.py', '3'))
report("W1", s3 and (s3['N'], s3['orbits'], s3['sizes']) == (219, 71, pred(3)['sizes']),
       "column-staged census over F_3: 219 folds in 71 G-orbits, orbit sizes 23×1, 11×2, 16×3, 21×6", s3['line'] if s3 else "no result")
t3 = parse(run('q123s2.py', '3'))
report("W2", t3 and s3 and all(t3[k] == s3[k] for k in ('N', 'orbits', 'sizes', 'nonfilt')) and matches(t3, 3),
       "torus-reduced census with incremental checks (q123s2.py) reproduces W1 exactly, including the fibre histogram",
       f"{t3['line'] if t3 else 'no result'} | fibres {t3['fib'] if t3 else ''}")
run('orbsizes.py', '3')
out = run('stabcheck.py')
want = "[((1, 'stab dim 2'), 23), ((2, 'unipotent'), 11), ((3, 'torus'), 16), ((6, 'stab dim 0'), 21)]"
report("W3", want in out and "by F_3 size: {}" in out,
       "every one of the 71 F_3 orbits has a member lifting to an exact Q-fold whose stabilizer in Lie(G) has the type its F_3 orbit "
       "size requires (1 ↔ all of G, 2 = q−1 ↔ unipotent line, 3 = q ↔ torus line, 6 = q(q−1) ↔ trivial)",
       " | ".join(l for l in out.splitlines() if 'orbits' in l))
ks = {}
for q in (3, 5, 7, 11, 13):
    m = re.search(r'algebra folds of K over F_\d+: (\d+)', run('kcount.py', str(q))); ks[q] = int(m.group(1)) if m else None
report("W4", all(ks[q] == 24 * q * q + 46 * q + 30 for q in ks),
       "the base of the kernel fibration (algebra folds of K = ker Δ) has 24q² + 46q + 30 points at q = 3, 5, 7, 11, 13", str(ks))
s5 = parse(run('q123s2.py', '5'))
report("W5", matches(s5, 5),
       "F_5 census: 567 = 21·25 + 30 + 12 folds in 71 orbits (23×1, 11×4, 16×5, 21×20), q − 1 = 4 non-filtered folds, and fibres "
       "{0: 24q²+19q+18, 1: 11q+12, q: 16q−2, q²: 1, 4q²−3q: 1} — the user's 9.7-hour run found the same 567 / 71",
       f"{s5['line'] if s5 else 'no result'} | fibres {s5['fib'] if s5 else ''}")
report("W6", s3 and s5 and s3['nonfilt'] == 2 and s5['nonfilt'] == 4,
       "the non-length-filtered folds (β ↦ β + δβ type) number q − 1 at q = 3, 5 (over F_3 they form one orbit of size 2)",
       f"{s3['nonfilt'] if s3 else '?'}, {s5['nonfilt'] if s5 else '?'}")
sp = run('special.py', '5')
report("W7", "{5: 1, 20: 4}" in sp.split("G-orbits on the special fibre: sizes")[-1].split("(")[0] if "G-orbits" in sp else False,
       "the special fibre (over P_K = id on vertices, δ ↦ 0) has 4q² − 3q points: four free G-orbits and one orbit of size q",
       " | ".join(l for l in sp.splitlines() if 'special fibre' in l))
if SLOW >= 1:
    s7 = parse(run('q123s2.py', '7'))
    sp7 = run('special.py', '7')
    report("W8", matches(s7, 7) and "{7: 1, 42: 4}" in sp7,
           "F_7 census: 1083 folds in 71 orbits (23×1, 11×6, 16×7, 21×42), 6 non-filtered, fibre histogram as predicted; special fibre 4 free + 1 of size 7",
           f"{s7['line'] if s7 else 'no result'} | fibres {s7['fib'] if s7 else ''}")
if SLOW >= 2:
    s11 = parse(run('q123s2.py', '11'))
    report("W9", matches(s11, 11), "F_11 census: 2619 folds in 71 orbits (23×1, 11×10, 16×11, 21×110), fibre histogram as predicted",
           f"{s11['line'] if s11 else 'no result'} | fibres {s11['fib'] if s11 else ''}")
print(f"\n{sum(res)}/{len(res)} checks pass")

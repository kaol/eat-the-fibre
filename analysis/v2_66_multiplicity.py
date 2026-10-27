"""v2.66: compare the vertex step of two quivers with the same support (one arrow per parallel class selected):
the strong Gröbner basis over ℤ of the vertex ideal (leading coefficient and monomial of every element), its vdim and
number of points over ℚ, and the point sets (mode p=0).
Usage: python3 analysis/v2_66_multiplicity.py Q1 Q2 OUT_DIR"""
import sys, os, re, subprocess, collections
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, ROOT)
import pcensus as pc
Q1, Q2, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
def run(Q, mode):
    f = os.path.join(out, f'{Q}_{mode.replace("=", "")}.txt')
    if not os.path.exists(f) or 'DONE' not in open(f).read():
        nv, A = pc.QUIVERS[Q]; cls = collections.OrderedDict()
        for i, st in enumerate(A): cls.setdefault(st, i)
        with open(f, 'w') as fh:
            subprocess.run(['nice', sys.executable, 'analysis/v2_58_vertex_integral.py', Q, os.path.join(out, f'{Q}_tmp'), mode,
                            '--arrow', ','.join(map(str, cls.values()))], stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
    return open(f).read()
def pts(t):
    res = set()
    for bl in t.split("POINTS")[1].split("--")[:-1]:
        body = bl.split("\n", 1)[1]
        res.add(tuple(sorted((m.group(1), -int(m.group(2)) if m.group(2) else 0) for m in re.finditer(r"(k\d\d)([-+]\d+)?\s*[,\n]", body + "\n"))))
    return res
r = {}
for Q in (Q1, Q2):
    ti, t0 = run(Q, 'int'), run(Q, 'p=0')
    r[Q] = dict(gb=sorted(re.findall(r"^LC (\S+ LM \S+)", ti, re.M)), vdim=re.search(r"VDIM (\d+)", ti).group(1), pts=pts(t0))
same = {k: r[Q1][k] == r[Q2][k] for k in ('gb', 'vdim', 'pts')}
print(f"{Q1} vs {Q2}: GB leading terms equal {same['gb']} ({len(r[Q1]['gb'])}, {len(r[Q2]['gb'])}); vdim {r[Q1]['vdim']}, {r[Q2]['vdim']}; "
      f"points {len(r[Q1]['pts'])}, {len(r[Q2]['pts'])}, equal {same['pts']}")

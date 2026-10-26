"""v2.57: tally of the v2.54 criterion over all vertex blocks of a quiver, from pipeline logs.

A run of a block is a header line "QUIVER V = BLOCK: presolved to … N components" followed by N component lines
("component k: dim … finitely many orbits: True" / "FAILS …" / "dim 0 (points)"); complete iff all N follow. In branch
mode (BRANCH=1) the block is "QUIVER V = BLOCK: L leaves", then L such runs, then "all leaves done; every leaf exited
cleanly"; complete iff all L runs are complete and that last line is present.
Per block, the outcome is
- PASS: some complete run whose components all pass or are points;
- FAIL: some complete run with a FAILS line (and no complete passing run);
- UNRESOLVED: no complete run.
Usage: python3 analysis/v2_57_tally.py QUIVER CENSUS_FILE LOG… → table on stdout."""
import re, sys
quiver, census, logs = sys.argv[1], sys.argv[2], sys.argv[3:]
blocks = []
for line in open(census):
    m = re.match(r"^([-0-9 ]+?)\s*\|", line)
    if not m or not re.search(r"\d", m.group(1)): continue
    blocks.append(",".join(r.replace('-1', 'm') for r in m.group(1).split()))
runs = {}                      # block -> list of (complete, n_components, n_pass, n_points, fails, leaves, source)
def add(b, rec): runs.setdefault(b, []).append(rec)
hdr = re.compile(rf"^{quiver} V = ([0-9m,]+): (.*)$")
for fn in logs:
    L = open(fn).read().splitlines()
    i = 0
    def unit(i):
        """parse one run starting at header line i; returns (next i, complete, ncomp, npass, npts, fails)"""
        m = hdr.match(L[i]); rest = m.group(2)
        mc = re.search(r"(\d+) components", rest)
        if not mc: return i + 1, False, 0, 0, 0, []
        n = int(mc.group(1)); j = i + 1; got = 0; npass = npts = 0; fails = []
        while j < len(L) and got < n and L[j].startswith("  component"):
            got += 1
            if "(points)" in L[j]: npts += 1
            elif "finitely many orbits: True" in L[j]: npass += 1
            elif "FAILS" in L[j]: fails.append(L[j].strip())
            j += 1
        return j, got == n, n, npass, npts, fails
    while i < len(L):
        m = hdr.match(L[i])
        if not m: i += 1; continue
        b, rest = m.group(1), m.group(2)
        ml = re.match(r"(\d+) leaves", rest)
        if ml:
            nl = int(ml.group(1)); i += 1; ok = True; nc = npass = npts = 0; fails = []; seen = 0
            while i < len(L) and seen < nl:
                if L[i].startswith("leaf "): i += 1; continue
                m2 = hdr.match(L[i])
                if not m2 or m2.group(1) != b: break
                i, c, n, p, q, f = unit(i); seen += 1; ok &= c; nc += n; npass += p; npts += q; fails += f
            while i < len(L) and L[i].startswith("leaf "): i += 1
            clean = i < len(L) and "all leaves done; every leaf exited cleanly" in L[i] and hdr.match(L[i]).group(1) == b
            add(b, (ok and seen == nl and clean, nc, npass, npts, fails, nl, fn)); i += 1
            continue
        if "all leaves done" in rest: i += 1; continue
        i, c, n, p, q, f = unit(i)
        add(b, (c, n, p, q, f, 0, fn))
out = {'PASS': 0, 'FAIL': 0, 'UNRESOLVED': 0}; ncomp = 0
for b in blocks:
    rs = runs.get(b, [])
    good = [r for r in rs if r[0] and not r[4] and r[2] + r[3] == r[1]]
    bad = [r for r in rs if r[0] and r[4]]
    if good:
        r = good[0]; st = 'PASS'; ncomp += r[1]
        det = f"{r[1]} {'leaf primes' if r[5] else 'components'} ({r[2]} pass, {r[3]} points)" + (f", {r[5]} leaves" if r[5] else "")
    elif bad:
        r = bad[0]; st = 'FAIL'; det = "; ".join(r[4])
    else:
        st = 'UNRESOLVED'; det = f"{len(rs)} incomplete runs"; r = (0,) * 6 + ('-',)
    out[st] += 1
    print(f"{b:24s} {st:10s} {det}  [{r[6].split('/')[-1]}]")
print(f"blocks {len(blocks)}: " + ", ".join(f"{k} {v}" for k, v in out.items()) + f"; components/leaf primes checked {ncomp}")

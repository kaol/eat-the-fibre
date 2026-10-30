"""v2.77: compare a v2_53_block_census table across primes by vertex block, identifying blocks mod 2 (the printed
symmetric lifts write 1 as −1 at p = 2; vertex blocks are pairwise distinct mod 2, v2.72). Prints the per-prime totals
and every block whose local-torus class count is not the same at all primes. (v2.78: any number of prime columns.)
Usage: python3 analysis/v2_77_blockdiff.py TABLE.txt"""
import sys, re, collections
def parse_row(s):
    return tuple(tuple(int(x) for x in re.findall(r'-?\d', r)) for r in s.split())
primes = None; cnt = collections.defaultdict(lambda: collections.Counter()); lift = {}
for line in open(sys.argv[1]):
    if line.startswith('vertex'):
        primes = [int(p) for p in re.findall(r'F_(\d+)', line)]; continue
    if '|' not in line: continue
    parts = [p.strip() for p in line.split('|')]
    B = parse_row(parts[0]); key = tuple(tuple(x % 2 for x in r) for r in B)
    for p, cell in zip(primes, parts[1:1 + len(primes)]):
        c = int(cell.split()[2]); cnt[key][p] += c
        if c and (p != 2 or key not in lift): lift[key] = parts[0]
print("totals:", {p: sum(cnt[k][p] for k in cnt) for p in primes}, f"blocks {len(cnt)}")
for k in sorted(cnt, key=lambda k: tuple(cnt[k][p] for p in primes)):
    vals = [cnt[k][p] for p in primes]
    if len(set(vals)) > 1: print(f"  {lift.get(k, k)}: " + ", ".join(f"q={p} {v}" for p, v in zip(primes, vals)))

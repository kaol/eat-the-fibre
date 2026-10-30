"""v2.77: compare a v2_53_block_census table at q = 2 and q = 3 by vertex block, identifying blocks mod 2 (the printed
symmetric lifts write 1 as −1 at p = 2; vertex blocks are pairwise distinct mod 2, v2.72). Prints the blocks whose
local-torus class counts differ.
Usage: python3 analysis/v2_77_blockdiff.py TABLE.txt"""
import sys, re, collections
def parse_row(s):
    return tuple(tuple(int(x) for x in re.findall(r'-?\d', r)) for r in s.split())
c2 = collections.Counter(); c3 = collections.Counter(); lift3 = {}
for line in open(sys.argv[1]):
    if '|' not in line or line.startswith('vertex'): continue
    parts = [p.strip() for p in line.split('|')]
    B = parse_row(parts[0]); f2 = list(map(int, parts[1].split())); f3 = list(map(int, parts[2].split()))
    key = tuple(tuple(x % 2 for x in r) for r in B)
    c2[key] += f2[2]; c3[key] += f3[2]
    if f3[2]: lift3[key] = parts[0]
tot2, tot3 = sum(c2.values()), sum(c3.values())
print(f"totals: q=2 {tot2}, q=3 {tot3}; blocks {len(set(c2) | set(c3))}")
for k in sorted(set(c2) | set(c3), key=lambda k: c3[k] - c2[k]):
    if c2[k] != c3[k]: print(f"  {lift3.get(k, k)}: q=2 {c2[k]}, q=3 {c3[k]} ({c3[k] - c2[k]:+d})")

"""v2.66: one row per quiver from logs/v2_66/*/summary.json, plus the q = 3 local-torus count of the census.
Usage: python3 analysis/v2_66_table.py"""
import json, glob, os, re
rows = []
for f in sorted(glob.glob('logs/v2_66/*/summary.json')):
    s = json.load(open(f)); Q = s['quiver']
    bc = f'logs/v2_66/{Q}.blockcensus.txt'
    loc = re.search(r"\{3: (\d+)\}", open(bc).read()).group(1) if os.path.exists(bc) else '—'
    t = re.match(r"blocks (\d+): PASS (\d+), FAIL (\d+), UNRESOLVED (\d+); components/leaf primes checked (\d+)", s['pipeline']['tally'])
    I = s['integral']; V = s['vertex']; G = s['goodprimes']
    bp = json.load(open(f'logs/v2_66/{Q}/vertex_badprimes.json'))
    vgood = all(all(r.values()) for r in bp['primes'].values())
    vbad = ','.join(map(str, bp['primes'])) or '—'
    ok = (I['snf1'] and I['saturated'] and I['detB'] in ('1', '-1') and I['group_identities'] and V['zero_dim'] and vgood
          and V['radical'] and V['point_ideal_vdim'] == V['points'] and t and t.group(1) == t.group(2)
          and set(map(int, s['pivots'])) <= {1, -1} and not s['branch_consts'] and G['exceptional'] == '[]' and G['cross_check_unresolved'] == 0)
    rows.append((Q, V['points'], vbad, t.group(5) if t else '?', len(s['pipeline']['case_split_blocks']), G['claims'],
                 ','.join(map(str, sorted(map(int, G['candidates'])))) or '—', s['bound_log10'], loc, ok))
print("| quiver | vertex blocks | vertex step checked directly at | components | case-split blocks | specialisation claims | candidate primes (all good) | log₁₀ B | classes at q = 3 | all checks |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows: print("| " + " | ".join(map(str, r)) + " |")

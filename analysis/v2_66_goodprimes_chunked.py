"""v2.66: v2_58_goodprimes.py in batches of blocks (one process per batch, so memory is returned between batches;
a single run on the 567 strata of 1⇉2→3→4→5 was killed for memory). Appends to OUT_TXT, skipping blocks already
reported there; prints EXCEPTIONAL as the union over batches.
Usage: python3 analysis/v2_66_goodprimes_chunked.py CERT OUT_PREFIX OUT_TXT [--batch N] [goodprimes options…]"""
import sys, os, json, subprocess, re
cert, outp, txt = sys.argv[1:4]; rest = sys.argv[4:]
n = 40
if rest[:1] == ['--batch']: n = int(rest[1]); rest = rest[2:]
blocks = []
for l in open(cert):
    r = json.loads(l)
    if r['kind'] == 'stratum' and r['block'] not in blocks: blocks.append(r['block'])
old = open(txt).read() if os.path.exists(txt) else ''
done = {m.group(1) for m in re.finditer(r"^(\S+) leaf \d+:", old, re.M)}
# point strata never print a leaf line: treat a block as done when every non-point stratum of it is reported
haspts = {}
for l in open(cert):
    r = json.loads(l)
    if r['kind'] == 'stratum': haspts.setdefault(r['block'], []).append(bool(r['vs']))
todo = [b for b in blocks if b not in done and any(haspts[b])]
exc = set(re.findall(r"^EXCEPTIONAL \[(.*)\]", old, re.M)[0].replace(' ', '').split(',')) - {''} if 'EXCEPTIONAL' in old else set()
with open(txt, 'a') as fh:
    for i in range(0, len(todo), n):
        batch = set(todo[i:i + n])
        skip = ';'.join(b for b in blocks if b not in batch)
        r = subprocess.run([sys.executable, 'analysis/v2_58_goodprimes.py', cert, outp, '--skip', skip] + rest, capture_output=True, text=True)
        out = r.stdout + r.stderr
        assert r.returncode == 0 and 'EXCEPTIONAL' in out, out[-2000:]
        exc |= set(re.search(r"EXCEPTIONAL \[(.*)\]", out).group(1).replace(' ', '').split(',')) - {''}
        fh.write("\n".join(l for l in out.splitlines() if not l.startswith('EXCEPTIONAL')) + "\n"); fh.flush()
        print(f"batch {i // n}: {len(batch)} blocks", flush=True)
    fh.write(f"EXCEPTIONAL [{', '.join(sorted(exc))}]\n")
print("EXCEPTIONAL", sorted(exc))

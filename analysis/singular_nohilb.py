"""v2.55: a copy of Singular's primdec.lib without Hilbert-driven Gröbner bases.

Singular 4.3.2's minAssGTZ overflows ("? overflow at t^24") in hilb() on the larger fold ideals (two_pairs_short), and
then silently returns an empty list. The Hilbert series there only drives std() (std(I, hilb, w) is std(I) with a hint).
This writes primdec.lib into DIR with:
- std(X, @hilb, @w), std(X, keephilb, @w) and std(X, hil, w) replaced by std(X);
- the hilb(…, 1, w) assignments replaced by intvec(0);
- groebner( replaced by std(.
Pointing SINGULARPATH at DIR makes every `LIB "primdec.lib"` (also the ones inside other libraries) load the copy. The
results are the same ideals; only the Gröbner strategy differs. Checked against v2.54 on pairs_in_short (v2.55 doc).
Usage: import singular_nohilb; singular_nohilb.install() → sets os.environ['SINGULARPATH']."""
import os, re
SRC = '/usr/share/singular/LIB/primdec.lib'
DIR = os.path.join(os.path.expanduser('~'), '.cache', 'etf_singular')
def install():
    os.makedirs(DIR, exist_ok=True)
    dst = os.path.join(DIR, 'primdec.lib')
    if not os.path.exists(dst):
        s = open(SRC).read()
        s = re.sub(r'std\(([^;]*?),\s*@hilb\s*,\s*@w\)', r'std(\1)', s)
        s = re.sub(r'std\(([^;]*?),\s*keephilb\s*,\s*@w\)', r'std(\1)', s)
        s = re.sub(r'std\(([^;]*?),\s*hil\s*,\s*w\)', r'std(\1)', s)
        s = re.sub(r'hilb\((\w+),\s*1,\s*@?w\)', r'intvec(0)', s)
        s = re.sub(r'\bgroebner\(', 'std(', s)
        open(dst, 'w').write(s)
    os.environ['SINGULARPATH'] = DIR
    return DIR

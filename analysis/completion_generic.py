"""δ-D2 for folds, generically.

Fix x ≺ y = x P(y) and x ≻ y = b1 xy + b2 P(x)y + b3 xP(y) + b4 P(x)P(y). Decompose each variable as p + n
(P p = p, P n = 0) and treat words in the pieces as free. P of a two-letter word is the word (im·im), 0 (ker·ker),
and for mixed words 0 or an opaque element, depending on the fold class:
  fold   : mixed words opaque          -> Baxter point only
  left   : ker·im ⊆ ker                -> line 1 + point 0
  right  : im·ker ⊆ ker                -> Baxter point only (≺ = xP(y) breaks the symmetry)
  bimod  : ker·im, im·ker ⊆ ker        -> v2 delta-D2: lines 1, 2 + points 0, (0,1,-1,0)
  + kik0 : also ker·im·ker = 0         -> the two points thicken to lines
Sufficiency is proved in lean/EatTheFibre/Discordia/Completion.lean; this script gives necessity for the free case.
Usage: python3 completion_generic.py {fold|left|bimod} [kik0]
"""
# Generic dendriform-completion classification for folds, in the decomposition x = p + n.
# Elements: dict word(tuple of letters) -> coeff. Letters: ('p',i),('n',i), or opaque ('U',a,b) = P(a b) unknown.
import sympy as sp, itertools, sys
b=sp.symbols('b1:5'); K=sp.Symbol('kik')  # kik: whether ker·im·ker products vanish (handled by mode)
def add(*es):
    r={}
    for e in es:
        for w,c in e.items(): r[w]=sp.expand(r.get(w,0)+c)
    return {w:c for w,c in r.items() if c!=0}
def sc(a,e): return {w:sp.expand(a*c) for w,c in e.items()}
def mul(e,f):
    r={}
    for (w1,c1),(w2,c2) in itertools.product(e.items(),f.items()):
        w=w1+w2; r[w]=sp.expand(r.get(w,0)+c1*c2)
    return {w:c for w,c in r.items() if c!=0}
def typ(l): return 'i' if l[0] in ('p','Q') else ('k' if l[0] in ('n','R') else '?')
MODE=sys.argv[1]   # 'fold' | 'left' (ker·im ⊆ ker) | 'right' (im·ker ⊆ ker) | 'bimod' (both)
KIK0=len(sys.argv)>2   # impose ker·im·ker = 0
def P(e):
    r={}
    for w,c in e.items():
        if len(w)==1:
            out={w:1} if typ(w[0])=='i' else {}
        elif len(w)==2:
            t=(typ(w[0]),typ(w[1]))
            if t==('i','i'): out={w:1}
            elif t==('k','k'): out={}
            elif t==('k','i') and MODE in('left','bimod'): out={}
            elif t==('i','k') and MODE in('right','bimod'): out={}
            else: out={(('Q',w),):1}   # opaque image element P(w); its complement is opaque kernel element
        else: raise Exception(w)
        r=add(r,sc(c,out))
    return r
def comp(e): return add(e,sc(-1,P(e)))
def el(i): return {(('p',i),):1,(('n',i),):1}
X,Y,Z=el(1),el(2),el(3)
def prec(x,y): return mul(x,P(y))
def succ(x,y): return add(sc(b[0],mul(x,y)),sc(b[1],mul(P(x),y)),sc(b[2],mul(x,P(y))),sc(b[3],mul(P(x),P(y))))
def star(x,y): return add(prec(x,y),succ(x,y))
def kill(e):
    if not KIK0: return e
    def bad(w):
        ts=[typ(l) for l in w]
        return any(ts[i]=='k' and ts[i+1]=='i' and ts[i+2]=='k' for i in range(len(ts)-2))
    return {w:c for w,c in e.items() if not bad(w)}
eqs=set()
for lhs,rhs in [(prec(prec(X,Y),Z),prec(X,star(Y,Z))),(prec(succ(X,Y),Z),succ(X,prec(Y,Z))),(succ(X,succ(Y,Z)),succ(star(X,Y),Z))]:
    for c in kill(add(lhs,sc(-1,rhs))).values(): eqs.add(c)
sol=sp.solve(list(eqs),b,dict=True)
print(MODE, 'kik=0' if KIK0 else '', sol)

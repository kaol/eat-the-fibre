"""δ-D2 completion variety of a census fold, read from the generated Lean data.

Usage: python3 completion_census.py lean/EatTheFibre/Gen/two_pairs.lean thm4P [t]
Results (v2.31): thm4P, pairs_in thm3F1P..F3P -> Baxter point only; thm3F4P -> line 1 + point 0;
two_pairs thm2P (an endomorphism, radical cube zero) -> four lines.
"""
import sys, sympy as sp, itertools
src=open(sys.argv[1]).read()
def grab(name):
    i=src.index(f"def {name}"); j=src.index(":=",i)+2
    k=src.find("\n\n",j); body=src[j:k] if k>0 else src[j:]
    return eval(body.split("/--")[0].replace("![","[").strip())
t=sp.Symbol('t'); b=sp.symbols('b1:5')
c=grab("c"); P=grab(sys.argv[2]); d=len(c)
Pm=sp.Matrix(d,d,lambda i,j: sum(a*t**n for n,a in enumerate(P[i][j])))
tv=sp.Rational(sys.argv[3]) if len(sys.argv)>3 else t
Pm=Pm.subs(t,tv)
def mul(u,w): return sp.Matrix([sum(c[i][j][m]*u[i]*w[j] for i in range(d) for j in range(d) if c[i][j][m]) for m in range(d)])
E=[sp.Matrix([1 if j==i else 0 for j in range(d)]) for i in range(d)]
pr=lambda x,y: mul(x,Pm*y)
su=lambda x,y: b[0]*mul(x,y)+b[1]*mul(Pm*x,y)+b[2]*mul(x,Pm*y)+b[3]*mul(Pm*x,Pm*y)
st=lambda x,y: pr(x,y)+su(x,y)
eqs=set()
for x,y,z in itertools.product(E,repeat=3):
    for e in list(pr(pr(x,y),z)-pr(x,st(y,z)))+list(pr(su(x,y),z)-su(x,pr(y,z)))+list(su(x,su(y,z))-su(st(x,y),z)):
        e=sp.expand(e)
        if e!=0: eqs.add(e)
print(len(eqs))
G=sp.groebner(list(eqs),*b,order='lex')
print(G.exprs)
print(sp.solve(list(G.exprs),b,dict=True))

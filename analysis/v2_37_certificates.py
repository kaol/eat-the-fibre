"""v2.37: integer certificates for the distinctness (rigidity) of the 19 families of 1⇉2⇉3 under G × T².

For each family P_t = Inn(1 + t·x) P₀ Inn(1 + t·x)⁻¹ (v2.36): every nonzero column of P_t lies in
V = ⟨e1, e2, e3, δa, δb, δaδb⟩ (asserted), on which any h ∈ G × T² acts by S_{l,m} = diag(1, 1, 1, l, m, lm).
If h P_t = P_t' h, then for each nonzero column j, with w_j = h e_j:  S_{l,m}(P_t e_j) = P_t' w_j, and P_t' w_j only
sees the coordinates of w_j at nonzero columns (the other columns of P_t' vanish). These are bilinear equations in
(w, l, m, t, t'). We search cofactors c_i (polynomials of degree ≤ D in w, l, m, t, t') with Σ c_i·eq_i = t' − t by a
sparse linear solve modulo a 31-bit prime, rational reconstruction, and an exact sympy check over ℚ[w, l, m, t, t'].
The certificates feed lean/gen/gen_rigid.py.
Usage: python3 analysis/v2_37_certificates.py NF.pkl REPS_two_pairs_7.pkl OUT.pkl [family ...]
"""
import sys, os, pickle, itertools, numpy as np, sympy as sp
from sympy.polys.monomials import itermonomials
nf, r7f, out = sys.argv[1:4]; which=[int(a) for a in sys.argv[4:]] or list(range(1,20))
HERE=os.path.dirname(os.path.abspath(__file__))
sys.argv=['x', r7f]; __file__=os.path.join(HERE,'family_normal_forms.py')
src=open(__file__).read()
exec(src.split("# ---- H = G × T² as matrices")[0]); exec(src[src.index("# ---- inner curves"):src.index("def is_fold_family")])
fams=pickle.load(open(nf,'rb')); t,tp,l,m=sp.symbols("t tp l m")
famlist=[(tag,i,xn) for tag in sorted(fams) for locs,(i,xn,_) in sorted(fams[tag].items())]
Vb=[[1 if q==idx[('v',0)] else 0 for q in range(d)],[1 if q==idx[('v',1)] else 0 for q in range(d)],[1 if q==idx[('v',2)] else 0 for q in range(d)],
    [{v(0):1,v(1):-1}.get(q,0) for q in range(d)],[{v(2):1,v(3):-1}.get(q,0) for q in range(d)],
    [{v(0,2):1,v(0,3):-1,v(1,2):-1,v(1,3):1}.get(q,0) for q in range(d)]]
Vm=sp.Matrix(Vb).T; scal=[1,1,1,l,m,l*m]
PR=2147483629
def modsolve(A,b):
    """solve A x = b mod PR (A: list of dict rows col->int), return dict or None"""
    n=max((c for row in A for c in row),default=-1)+1
    rows=[(dict(row),bb%PR) for row,bb in zip(A,b)]
    piv={}  # col -> (row dict, rhs)
    for row,bb in rows:
        row={c:v%PR for c,v in row.items() if v%PR}
        for c in sorted(list(row)):
            if c in piv and c in row:
                f=row[c]; pr,pb=piv[c]
                for cc,vv in pr.items(): row[cc]=(row.get(cc,0)-f*vv)%PR
                bb=(bb-f*pb)%PR; row={cc:vv for cc,vv in row.items() if vv}
        if not row:
            if bb: return None
            continue
        c=min(row); inv=pow(row[c],PR-2,PR)
        row={cc:vv*inv%PR for cc,vv in row.items()}; bb=bb*inv%PR
        # eliminate c from existing pivots
        for c2,(pr,pb) in list(piv.items()):
            if c in pr:
                f=pr[c]; nr=dict(pr)
                for cc,vv in row.items(): nr[cc]=(nr.get(cc,0)-f*vv)%PR
                nr={cc:vv for cc,vv in nr.items() if vv}; piv[c2]=(nr,(pb-f*bb)%PR)
        piv[c]=(row,bb)
    x={}
    for c,(row,bb) in piv.items(): x[c]=bb   # free vars = 0 (rows are fully reduced w.r.t. other pivots)
    return x
def ratrec(a):
    a%=PR; r0,r1=PR,a; s0,s1=0,1
    while r1*r1*2>PR: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    return sp.Rational(r1,s1) if s1 else None
results={}
for k in which:
    tag,i,xn=famlist[k-1]
    Pz=((r['greps'][i].astype(np.int64)%7+3)%7)-3; Inn,Innv=INN[xn]
    P=(Inn*sp.Matrix(Pz.tolist())*Innv).applyfunc(sp.expand).subs(s_,t)
    nz=[j for j in range(d) if any(P[q,j]!=0 for q in range(d))]
    eqs=[]; xs=[]; cols={}
    for j in nz:
        c=sp.symbols('c0:6'); s=sp.solve(list(Vm*sp.Matrix(c)-P[:,j]),c,dict=True)
        assert s, f"family {k}: column {names[j]} not in V"
        coords=[sp.expand(s[0][ci]) for ci in c]; cols[j]=coords
        Scol=Vm*sp.Matrix([sc*co for sc,co in zip(scal,coords)])
        w=[sp.Symbol(f'w{j}_{ii}') for ii in nz]; xs+=w
        rhs=sum((wi*P[:,ii].subs(t,tp) for wi,ii in zip(w,nz)),sp.zeros(d,1))
        for q in range(d):
            e=sp.expand(Scol[q]-rhs[q])
            if e!=0: eqs.append((j,q,e))
    gens=xs+[l,m,t,tp]
    found=None
    for D in (1,2):
        mons=sorted(itermonomials(gens,D),key=sp.default_sort_key)
        unk=[(ei,mo) for ei in range(len(eqs)) for mo in mons]
        # Σ c_{ei,mo} * mo * e_ei = tp - t   (coefficientwise)
        coeff_rows={}
        for ui,(ei,mo) in enumerate(unk):
            for mono,cf in sp.Poly(sp.expand(mo*eqs[ei][2]),*gens).terms():
                coeff_rows.setdefault(mono,{})[ui]=int(cf)
        target=sp.Poly(tp-t,*gens).as_dict()
        keys=list(coeff_rows)
        for mono in target:
            if mono not in coeff_rows: coeff_rows[mono]={}; keys.append(mono)
        A=[coeff_rows[mo] for mo in keys]; b=[int(target.get(mo,0)) for mo in keys]
        sol=modsolve(A,b)
        if sol is None: continue
        cof=[0]*len(eqs)
        for ui,val in sol.items():
            rv=ratrec(val); ei,mo=unk[ui]; cof[ei]+=rv*mo
        chk=sp.expand(sum(c*e for c,(_,_,e) in zip(cof,eqs))-(tp-t))
        if chk==0: found=(D,cof); break
        print('family',k,'D',D,'reconstruction failed verification')
    print(f'family {k}: x {xn}, nz {[names[j] for j in nz]}, eqs {len(eqs)}, certificate degree {found[0] if found else None}',flush=True)
    results[k]=dict(tag=tag,rep=i,x=xn,nz=nz,cols={j:[str(c) for c in cols[j]] for j in nz},
                    eqs=[(j,q,str(e)) for j,q,e in eqs],cof=[str(sp.expand(c)) for c in found[1]] if found else None)
pickle.dump(results,open(out,'wb'))
print('denominators:', sorted({sp.fraction(sp.together(sp.sympify(c)))[1] for v_ in results.values() for c in v_['cof']}, key=str))

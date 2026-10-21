import sys, pickle, collections, numpy as np, itertools
files=sys.argv[1:]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'family_signatures.py')
exec(open(__file__).read().split("\ndef main():")[0])
import os; FR={'__file__':os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'fibre_rigidity.py')}
exec(open(FR['__file__']).read().split("if __name__ == '__main__':")[0].replace("_files = sys.argv[1:]; sys.argv = ['x']; _cwd","_cwd"), FR)
rank,nullspace,jacobian=FR['rank'],FR['nullspace'],FR['jacobian']
out={}
for f in files:
    r=pickle.load(open(f,'rb')); name,p=r['name'],r['p']
    nv,arrows=QUIVERS[name]; basis,M,DL=build(nv,arrows); M=np.array(M,np.int64); d=len(basis)
    C=FR['setup'](name,p); idx={b:i for i,b in enumerate(basis)}
    v=lambda *path: idx[('p',path)]
    da=np.zeros(d,np.int64); da[v(0)]=1; da[v(1)]=-1
    db=np.zeros(d,np.int64); db[v(2)]=1; db[v(3)]=-1
    mulv=lambda x,y: np.einsum('i,j,ijm->m',x,y,M)
    dadb=mulv(da,db)%p
    def ad(x):  # matrix of y ↦ xy − yx
        L=np.einsum('i,ijm->mj',x,M); R=np.einsum('j,ijm->mi',x,M)
        return (L-R)%p
    cands={'ad(δaδb)':ad(dadb),'ad(δa)':ad(da%p),'ad(δb)':ad(db%p),'E_a':C['E'][0],'E_b':C['E'][1]}
    names=list(cands)
    loc=r['loc']; seen=set(); grp=collections.defaultdict(collections.Counter)
    for i,P in enumerate(r['greps']):
        if loc[i] in seen: continue
        seen.add(loc[i]); P=P.astype(np.int64)%p
        s=signature(P,M,p)
        J=jacobian(P,C,p); br=lambda D:((D@P-P@D)%p).ravel()
        # S
        cond=(np.einsum('ml,cmab->clab',P,C['DLc'])-np.einsum('ai,clij,bj->clab',P,C['DLc'],P))%p
        Wl=nullspace(cond.reshape(C['ncl'],-1).T%p,p)
        Sd=[br(D) for D in C['DerD']]+[br(np.einsum('c,cij->ij',w,C['E'])%p) for w in Wl]
        Sm=np.array(Sd) if Sd else np.zeros((0,d*d),np.int64); rS=rank(Sm,p) if len(Sm) else 0
        # motif subspace V: combos c with [c,P] ∈ T
        A=np.array([(J@br(cands[n]))%p for n in names]).T
        Wv=nullspace(A,p)
        vecs=[br(sum(int(c)*cands[n] for c,n in zip(w,names))%p) for w in Wv]
        ext=rank(np.concatenate([Sm,np.array(vecs)]),p)-rS if len(vecs) else 0
        # which single motifs (or pairs) realise excess
        used=tuple(n for k,n in enumerate(names) if any(w[k] for w in Wv)) if ext else ()
        grp[s][(ext,used)]+=1
    out[p]=grp
p1,p2=sorted(out)
for s in sorted(set(out[p1])|set(out[p2]),key=lambda s:(s[0],str(s))):
    a,b=out[p1].get(s,collections.Counter()),out[p2].get(s,collections.Counter())
    if sum(a.values())!=sum(b.values()):
        print(f"{s[0]:5} im {s[2]} ker {s[3]}:")
        for k in sorted(set(a)|set(b),key=str):
            if a[k]!=b[k] or k[0]: print(f"      motif excess {k[0]} via {k[1]}: {a[k]}→{b[k]}")

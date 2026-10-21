import sys, pickle, collections, numpy as np
files=sys.argv[1:]; sys.argv=['x']
import os; FR={'__file__':os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'fibre_rigidity.py')}
exec(open(FR['__file__']).read().split("if __name__ == '__main__':")[0].replace("_files = sys.argv[1:]; sys.argv = ['x']; _cwd","_cwd"), FR)
rank,nullspace,jacobian,residual_series,setup=FR['rank'],FR['nullspace'],FR['jacobian'],FR['residual_series'],FR['setup']
QUIVERS,build=FR['QUIVERS'],FR['build']
for f in files:
    r=pickle.load(open(f,'rb')); name,p=r['name'],r['p']
    nv,arrows=QUIVERS[name]; basis,M,DL=build(nv,arrows); M=np.array(M,np.int64); d=len(basis)
    C=setup(name,p); idx={b:i for i,b in enumerate(basis)}
    cls=[[i for i,a in enumerate(arrows) if a==st] for st in sorted(set(arrows)) if arrows.count(st)==2]
    (a1,a2),(b1,b2)=cls
    w=np.zeros(d,np.int64); w[idx[('p',(b1,))]]=1; w[idx[('p',(b2,))]]=-1; w[idx[('p',(a1,))]]-=1; w[idx[('p',(a2,))]]+=1
    w%=p
    D=(np.einsum('i,ijm->mj',w,M)-np.einsum('j,ijm->mi',w,M))%p
    loc=r['loc']; seen=set(); cnt=collections.Counter()
    inv2=pow(2,p-2,p)
    for i,P in enumerate(r['greps']):
        if loc[i] in seen: continue
        seen.add(loc[i]); P=P.astype(np.int64)%p
        J=jacobian(P,C,p); br=lambda X:((X@P-P@X)%p)
        V=br(D)
        if not V.any(): cnt['stabilises P']+=1; continue
        if (J@V.ravel()%p).any(): cnt['not tangent']+=1; continue
        cond=(np.einsum('ml,cmab->clab',P,C['DLc'])-np.einsum('ai,clij,bj->clab',P,C['DLc'],P))%p
        Wl=nullspace(cond.reshape(C['ncl'],-1).T%p,p)
        Sd=[br(X).ravel() for X in C['DerD']]+[br(np.einsum('c,cij->ij',x,C['E'])%p).ravel() for x in Wl]
        Sm=np.array(Sd) if Sd else np.zeros((0,d*d),np.int64)
        if len(Sm) and rank(np.concatenate([Sm,V.ravel()[None]]),p)==rank(Sm,p): cnt['census symmetry']+=1; continue
        imgJ=np.array([(J@br(D2).ravel())%p for D2 in C['Der']]).T
        W0=(D@V-V@D)%p*inv2%p
        r2=residual_series(np.array([P,V,W0]),C,p)
        ok=rank(np.concatenate([imgJ,((-r2)%p)[:,None]],1),p)==rank(imgJ,p)
        cnt['excess, extends to 2nd order' if ok else 'excess, obstructed at 2nd order']+=1
    print(name,p,dict(sorted(cnt.items())))

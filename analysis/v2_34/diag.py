import sys, pickle, numpy as np, collections
f=sys.argv[1]; sys.argv=['x',f]
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
src=open(__file__).read().split("def main")[0]; exec(src)
r=pickle.load(open(f,'rb')); name,p=r['name'],r['p']; C=setup(name,p); d=C['d']
cnt=collections.Counter()
for i,P in enumerate(r['greps']):
    P=P.astype(np.int64)%p; a=at(P,C,p)
    if not a['e']: continue
    J=jacobian(P,C,p); br=lambda D: ((D@P-P@D)%p).ravel()
    # torus directions tangent to Folds at first order
    Etan=[w for w in np.eye(C['ncl'],dtype=np.int64)]
    # all w with [E_w,P] in T: nullspace over w
    A=np.array([J@br(np.einsum('c,cij->ij',w,C['E'])%p)%p for w in np.eye(C['ncl'],dtype=np.int64)]).T
    Wfirst=nullspace(A%p,p)
    cond=(np.einsum('ml,cmab->clab',P,C['DLc'])-np.einsum('ai,clij,bj->clab',P,C['DLc'],P))%p
    W=nullspace(cond.reshape(C['ncl'],-1).T%p,p)
    S=[br(D) for D in C['DerD']]+[br(np.einsum('c,cij->ij',w,C['E'])%p) for w in Wfirst]
    dimS2=rank(np.array(S),p) if S else 0
    cnt[(a['e'], len(Wfirst)-len(W), a['F']-dimS2)]+=1
print(name,p,'(e, extra first-order torus weights, e after adding them):',dict(cnt))

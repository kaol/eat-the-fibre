import sys, pickle, numpy as np, collections, itertools
q, rfile = sys.argv[1:3]; extra=sys.argv[3:]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
r=pickle.load(open(rfile,'rb')); name,p=r['name'],r['p']; C=setup(name,p); d=C['d']
# sanity: Der_Δ directions extend at second order at every rep (they are genuine curves)
bad=0
for P in r['greps'][:40]:
    P=P.astype(np.int64)%p; J=jacobian(P,C,p); br=lambda D:(D@P-P@D)%p
    imgJ=np.array([(J@br(D2).ravel())%p for D2 in C['Der']]).T; rimg=rank(imgJ,p)
    for D in C['DerD']:
        V=br(D); W0=(D@V-V@D)%p*pow(2,p-2,p)%p
        r2=residual_series(np.array([P,V,W0]),C,p)
        if rank(np.concatenate([imgJ,((-r2)%p)[:,None]],1),p)!=rimg: bad+=1
print(q,'sanity: Der_Δ directions failing the second-order test:',bad)
res=collections.Counter()
for i,P in enumerate(r['greps']):
    e,ext=second_order(P,C,p)
    if e: res[(e,ext)]+=1
print(q,'(e, extendable excess classes) over reps with e>0:',dict(res))

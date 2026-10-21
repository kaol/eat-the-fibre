import sys, pickle, numpy as np, collections
rfile, ffile, bfile, out = sys.argv[1:5]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
r=pickle.load(open(rfile,'rb')); name,p=r['name'],r['p']; C=setup(name,p)
folds=np.array(pickle.load(open(ffile,'rb'))['folds'],np.int64)%p
B=pickle.load(open(bfile,'rb')); comp=B['comp0']
key={f.astype(np.int8).tobytes():i for i,f in enumerate(folds)}
fib=collections.defaultdict(list)
for i,P in enumerate(r['greps']):
    j=key[(P.astype(np.int64)%p).astype(np.int8).tobytes()]
    a=at(P,C,p); e,ext=second_order(P,C,p); nproj=sum(p**k for k in range(e))
    eff=0 if ext==0 else (e if ext==nproj else 1)
    fib[comp[j]].append((i,a['S'],a['S']+eff))
sel=[]
for b,L in fib.items():
    mS=max(x[1] for x in L)
    if max(x[2] for x in L)>mS: sel.append(max(L,key=lambda x:x[2])[0])
pickle.dump(sel,open(out+'.sel','wb')); print(len(sel),'excess classes; representatives',sel)

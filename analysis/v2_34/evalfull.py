import sys, pickle, numpy as np, collections
rfile, ffile, bfile, resf = sys.argv[1:5]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
r=pickle.load(open(rfile,'rb')); name,p=r['name'],r['p']; C=setup(name,p)
folds=np.array(pickle.load(open(ffile,'rb'))['folds'],np.int64)%p
B=pickle.load(open(bfile,'rb')); Gorb=B['Gorb']
key={f.astype(np.int8).tobytes():i for i,f in enumerate(folds)}
g2rep={}
for i,P in enumerate(r['greps']): g2rep[int(Gorb[key[(P.astype(np.int64)%p).astype(np.int8).tobytes()]])]=i
res=pickle.load(open(resf,'rb'))
cache={}
def val(i):
    if i not in cache:
        P=r['greps'][i]; a=at(P,C,p); e,ext=second_order(P,C,p); nproj=sum(p**k for k in range(e))
        eff=0 if ext==0 else (e if ext==nproj else 1); cache[i]=(a['S'],a['S']+eff)
    return cache[i]
done=set(); out=collections.Counter()
for gi,hits in res.items():
    beh=frozenset(hits)
    if beh in done: continue
    done.add(beh)
    vals=[val(g2rep[g]) for g in hits]
    mS=max(v[0] for v in vals); m2=max(v[1] for v in vals)
    out[f'second-order excess {m2-mS}']+=1
    print(f'behaviour of rep {gi}: {len(hits)} G-orbits, max S = {mS}, max second-order = {m2}')
print(dict(out))

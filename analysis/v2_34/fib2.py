import sys, pickle, numpy as np, collections, itertools
q, rfile, ffile, bfile = sys.argv[1:5]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
r=pickle.load(open(rfile,'rb')); name,p=r['name'],r['p']; C=setup(name,p)
folds=np.array(pickle.load(open(ffile,'rb'))['folds'],np.int64)%p
B=pickle.load(open(bfile,'rb')); comp=B['comp0']; Gorb=B['Gorb']
key={f.astype(np.int8).tobytes():i for i,f in enumerate(folds)}
fib=collections.defaultdict(list)
for i,P in enumerate(r['greps']):
    j=key[(P.astype(np.int64)%p).astype(np.int8).tobytes()]
    a=at(P,C,p); e,ext=second_order(P,C,p)
    nproj=sum(p**k for k in range(e))            # number of projective classes in F/S
    eff = 0 if ext==0 else (e if ext==nproj else 1)
    fib[comp[j]].append((a['S'], a['S']+e, a['S']+eff))
# width of each behaviour in G-orbits
width=collections.Counter()
seen=set()
for i in range(len(folds)):
    if (comp[i],Gorb[i]) not in seen: seen.add((comp[i],Gorb[i])); width[comp[i]]+=1
stat=collections.Counter()
for b,L in fib.items():
    mS=max(x[0] for x in L); m1=max(x[1] for x in L); m2=max(x[2] for x in L)
    stat[(f'first-order excess {m1-mS}', f'second-order excess {m2-mS}', f'width {width[b]}' if width[b]>4 else 'width ≤4')]+=1
print(q)
for k,v in sorted(stat.items()): print('   ',v,k)

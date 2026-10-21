"""Full Aut(kQ)-orbit of selected representatives on 1⇉2⇇3 + 1→3, including a_i ↦ a_i + Σ_j λ_ij c b_j."""
import sys, pickle, itertools, collections, numpy as np
rfile, ffile, bfile, out = sys.argv[1:5]; sys.argv=['x']
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fibre_rigidity.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
name,p='pairs_in_short',3
r=pickle.load(open(rfile,'rb')); C=setup(name,p); d=C['d']
nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); M=np.array(M,np.int64)
idxb={b:i for i,b in enumerate(basis)}
folds=np.array(pickle.load(open(ffile,'rb'))['folds'],np.int64)%p; N=len(folds)
B=pickle.load(open(bfile,'rb')); comp=B['comp0']; Gorb=B['Gorb']
rng=np.random.default_rng(1); W=rng.integers(1,2**62,size=d*d,dtype=np.int64)
hsh=lambda X:(X.reshape(len(X),-1)*W).sum(1); H=hsh(folds); order=np.argsort(H); Hs=H[order]
def lookup(X):
    h=hsh(X); pos=np.searchsorted(Hs,h).clip(0,N-1); idx=order[pos]
    return np.where((Hs[pos]==h)&(folds[idx]==X).all((1,2)),idx,-1)
classes=collections.defaultdict(list)
for ai,st in enumerate(arrows): classes[st].append(ai)
cl=list(classes.values())
def auto(img_ar):
    """algebra endomorphism fixing vertices with arrow x ↦ img_ar[x] (vector); paths multiplicatively"""
    mulv=lambda x,y: np.einsum('i,j,ijm->m',x,y,M)%p
    Phi=np.zeros((d,d),np.int64)
    for i,b in enumerate(basis):
        if b[0]=='v': Phi[i,i]=1; continue
        v=img_ar[b[1][0]]
        for x in b[1][1:]: v=mulv(v,img_ar[x])
        Phi[:,i]=v
    return Phi
def inv(A):
    import sympy as sp
    return np.array(sp.Matrix(A.tolist()).inv_mod(p).tolist(),np.int64)
E=lambda i: np.eye(d,dtype=np.int64)[i]
arrowidx={ai: idxb[('p',(ai,))] for ai in range(len(arrows))}
# coset reps of G in Sub (one per column-sum vector per class)
def gl(m):
    for vals in itertools.product(range(p),repeat=m*m):
        A=np.array(vals,np.int64).reshape(m,m)
        if round(np.linalg.det(A))%p: yield A
cos=[]
for c in cl:
    seen={}
    for A in gl(len(c)): seen.setdefault(tuple(A.sum(0)%p),A)
    cos.append(list(seen.values()))
# parallel paths to arrows: for arrow x (s→t), paths of length ≥ 2 from s to t
par={x:[i for i,b in enumerate(basis) if b[0]=='p' and len(b[1])>=2 and arrows[b[1][0]][0]==arrows[x][0] and arrows[b[1][-1]][1]==arrows[x][1]] for x in range(len(arrows))}
extra=[(x,i) for x in range(len(arrows)) for i in par[x]]
print('higher substitutions:',[(x,basis[i]) for x,i in extra])
Ss=[]
for Ms in itertools.product(*cos):
    for lam in itertools.product(range(p),repeat=len(extra)):
        img={}
        for Mc,c in zip(Ms,cl):
            for jx,x in enumerate(c):
                v=np.zeros(d,np.int64)
                for iy,y in enumerate(c): v[arrowidx[y]]=Mc[iy,jx]%p
                img[x]=v
        for (x,i),l in zip(extra,lam): img[x]=img[x].copy(); img[x][i]=(img[x][i]+l)%p
        Ss.append(auto(img))
Sm=np.array(Ss,np.float64); Si=np.array([inv(s) for s in Ss],np.float64)
print('substitution coset reps:',len(Ss))
vert=[i for i,b in enumerate(basis) if b[0]=='v']; rad=[i for i,b in enumerate(basis) if b[0]=='p']
one=np.zeros(d,np.int64); one[vert]=1
mulb=lambda x,y: np.einsum('...i,...j,ijm->...m',x,y,M)%p
R=np.array(list(itertools.product(range(p),repeat=len(rad))),np.int64); rr=np.zeros((len(R),d),np.int64); rr[:,rad]=R
u=(one+rr)%p; pw=np.broadcast_to(one,rr.shape).copy(); uinv=pw.copy()
for _ in range(d): pw=mulb(pw,(-rr)%p); uinv=(uinv+pw)%p
def innmat(a,b):
    La=np.einsum('ni,ijm->nmj',a,M)%p; Rb=np.einsum('nj,ijm->nmi',b,M)%p
    return np.einsum('nmk,nkj->nmj',Rb,La)%p
Inn=innmat(u,uinv).astype(np.float64); Innv=innmat(uinv,u).astype(np.float64)
key={f.astype(np.int8).tobytes():i for i,f in enumerate(folds)}
targets=[int(x) for x in open(sys.argv[0].replace('fullorbit.py','targets.txt')).read().split()] if False else None
reps_idx=[key[(P.astype(np.int64)%p).astype(np.int8).tobytes()] for P in r['greps']]
sel=pickle.load(open(out+'.sel','rb'))
res={}
CH=max(1,(1<<27)//(len(Inn)*d*d*8))
for gi in sel:
    j=reps_idx[gi]; P=folds[j].astype(np.float64); hits=set()
    for c0 in range(0,len(Ss),CH):
        X=np.mod(Sm[c0:c0+CH]@P@Si[c0:c0+CH],p)
        Q=np.mod(Inn[None]@X[:,None]@Innv[None],p).astype(np.int64).reshape(-1,d,d)
        h=lookup(Q); h=h[h>=0]
        hits|=set(int(Gorb[x]) for x in np.unique(h))
    res[gi]=sorted(hits); print(gi,'reaches',len(hits),'G-orbits; classes',sorted({int(comp[x]) for x in range(N) if Gorb[x] in hits}) if len(hits)<50 else '…',flush=True)
pickle.dump(res,open(out,'wb'))

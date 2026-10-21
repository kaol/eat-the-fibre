import sys, pickle, collections, numpy as np
files=sys.argv[1:]; sys.argv=['x']
import os
import os; __file__=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'family_signatures.py')
src=open(__file__).read().split("\ndef main():")[0]; exec(src)
import os; FR={'__file__':os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'fibre_rigidity.py')}; exec(open('/home/claude/eat-the-fibre/analysis/fibre_rigidity.py').read().split("if __name__ == '__main__':")[0].replace("_files = sys.argv[1:]; sys.argv = ['x']; _cwd","_cwd"), FR)
out={}
for f in files:
    r=pickle.load(open(f,'rb')); name,p=r['name'],r['p']
    nv,arrows=QUIVERS[name]; basis,M,DL=build(nv,arrows); M=np.array(M,np.int64)
    C=FR['setup'](name,p); loc=r['loc']
    grp=collections.defaultdict(lambda:[0,0]); seen=set()
    for i,P in enumerate(r['greps']):
        if loc[i] in seen: continue
        seen.add(loc[i])
        s=signature(P,M,p); e,ext=FR['second_order'](P,C,p)
        g=grp[s]; g[0]+=1; g[1]+= 1 if ext>0 else 0
    out[p]=grp
p1,p2=sorted(out)
for s in sorted(set(out[p1])|set(out[p2]),key=lambda s:(s[0],str(s))):
    a=out[p1].get(s,[0,0]); b=out[p2].get(s,[0,0])
    if a[0]!=b[0]: print(f"{s[0]:5} im {s[2]} ker {s[3]}: local orbits {a[0]}→{b[0]}, with extendable excess {a[1]}→{b[1]}")
tot=[sum(v[1] for v in out[p].values()) for p in (p1,p2)]
print('local orbits with extendable excess, all signatures:',tot)

#include <stdio.h>
#include <stdlib.h>
/* forests: par[v] in 0..n-1 (point) or n (R1) or n+1 (R2); acyclic.
   node v = block {v} ∪ children(v); side = root side if depth even (depth = #edges to root, >=1) ... match python:
   depth_root: d counts steps until reaching root; side = rs if d%2==0 else 3-rs.
   S_v = union over x in subtree(v) with side[x]==1 of {x} ∪ children(x).
   compat: for each i: S_i empty or S_{i+1} empty or (S_i={k}, S_{i+1}={k+1}). */
int n; int par[16];
long long total=0, good=0;
int main(int argc,char**argv){
  n=atoi(argv[1]); int B=n+2;
  for(int i=0;i<n;i++)par[i]=0;
  long long cnt[64]={0};
  while(1){
    /* acyclicity + depth/root */
    int depth[16],root[16],ok=1;
    for(int v=0;v<n;v++){ int x=v,d=0; while(x<n && d<=n){x=par[x];d++;} if(d>n){ok=0;break;} depth[v]=d; root[v]=x; }
    if(ok){
      total++;
      int side[16]; unsigned ch[16]; for(int v=0;v<n;v++){ch[v]=0;}
      for(int v=0;v<n;v++){ int rs=(root[v]==n)?1:2; side[v]=(depth[v]%2==0)?rs:3-rs; if(par[v]<n) ch[par[v]]|=1u<<v; }
      /* S_v via post-order: compute subtree masks: process by decreasing depth */
      unsigned S[16]; int order[16]; for(int i=0;i<n;i++)order[i]=i;
      for(int i=0;i<n;i++)for(int j=i+1;j<n;j++) if(depth[order[j]]>depth[order[i]]){int t=order[i];order[i]=order[j];order[j]=t;}
      for(int i=0;i<n;i++){int v=order[i]; unsigned s = (side[v]==1)? ((1u<<v)|ch[v]) : 0u;
        for(int c=0;c<n;c++) if(ch[v]>>c&1) s|=S[c]; S[v]=s;}
      int c=1;
      for(int i=0;i+1<n && c;i++){ unsigned a=S[i],b=S[i+1]; if(!a||!b) continue;
        if(__builtin_popcount(a)==1 && __builtin_popcount(b)==1 && b==(a<<1)) continue; c=0; }
      good+=c;
    }
    int k=0; while(k<n){ par[k]++; if(par[k]<B) break; par[k]=0; k++; } if(k==n) break;
  }
  printf("n=%d forests=%lld compatible=%lld\n",n,total,good);
}

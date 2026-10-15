import sys, os, json, pathlib, numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import *
from wr_model import haversine, Ship
from build_mesh import *
G=get_grid()
M=Mesh(G)
STATS={}
def make(key, prefix):
    src=pd.read_csv(DATA/f'{key}_nos_consultados.csv')
    q=src[src['coluna']==src['coluna'].round()].copy().reset_index(drop=True)      # so nos com clima consultado
    q['coluna']=q['coluna'].round().astype(int)
    orig=q[(q.coluna==q.coluna.min())&(q.offset_nm==0)].iloc[0]; dest=q[(q.coluna==9)&(q.offset_nm==0)].iloc[0]
    ports=[(orig.Latitude,orig.Longitude),(dest.Latitude,dest.Longitude)]
    # 1) folga dos nos
    q['clr']=G.clearance(q.Latitude.values,q.Longitude.values)
    q['near_port']=[any(haversine(a,b,pl,po)<=PORT_MARGIN for pl,po in ports) for a,b in zip(q.Latitude,q.Longitude)]
    q['land']=globe.is_land(q.Latitude.values,q.Longitude.values) & ~np.array([any(haversine(a,b,pl,po)<=PORT_TOL for pl,po in ports) for a,b in zip(q.Latitude,q.Longitude)])
    q['in_mesh']=((q.clr>=D_MIN)|q.near_port)&(~q.land)
    q.loc[q.id.isin([orig.id,dest.id]),'in_mesh']=True
    nid=int(src.id.max())+1
    nodes=q.copy(); newn=[]; arcs=[]
    st=dict(nodes_q=len(q), nodes_dropped=int((~q.in_mesh).sum()), pairs=0, chords_ok=0, planned=0, dropped=0, waypoints=0, ratios=[])
    cols=sorted(q.coluna.unique())
    for c0,c1 in zip(cols[:-1],cols[1:]):
        A=q[(q.coluna==c0)&q.in_mesh]; B=q[(q.coluna==c1)&q.in_mesh]
        for a in A.itertuples():
            for b in B.itertuples():
                st['pairs']+=1
                ok,mc=M.seg_ok(a.Latitude,a.Longitude,b.Latitude,b.Longitude,ports)
                if ok: arcs.append((int(a.id),int(b.id))); st['chords_ok']+=1; continue
                Q=M.plan(a.Latitude,a.Longitude,b.Latitude,b.Longitude,ports)
                if Q is None: st['dropped']+=1; continue
                L=sum(float(haversine(Q[i,0],Q[i,1],Q[i+1,0],Q[i+1,1])) for i in range(len(Q)-1))
                ch=float(haversine(a.Latitude,a.Longitude,b.Latitude,b.Longitude))
                if L>MAX_RATIO*ch: st['dropped']+=1; continue
                st['planned']+=1; st['ratios'].append(L/ch)
                cum=np.r_[0,np.cumsum([haversine(Q[i,0],Q[i,1],Q[i+1,0],Q[i+1,1]) for i in range(len(Q)-1)])]
                prev=int(a.id)
                for k in range(1,len(Q)-1):
                    row={c:np.nan for c in nodes.columns}
                    row.update(id=nid,Latitude=Q[k,0],Longitude=Q[k,1],coluna=a.coluna+(b.coluna-a.coluna)*cum[k]/cum[-1],offset_nm=0,in_mesh=True,chain_from=int(a.id),chain_to=int(b.id),chain_frac=float(cum[k]/cum[-1]))
                    newn.append(row); arcs.append((prev,nid)); prev=nid; nid+=1; st['waypoints']+=1
                arcs.append((prev,int(b.id)))
    nodes=pd.concat([nodes,pd.DataFrame(newn)],ignore_index=True)
    nodes.drop(columns=['clr','near_port','land','in_mesh'],errors='ignore').to_csv(MESHDIR/f'{key}_nos.csv',index=False)
    pd.DataFrame(arcs,columns=['origem','destino']).drop_duplicates().to_csv(MESHDIR/f'{key}_arcos.csv',index=False)
    st['arcs']=len(set(arcs)); st['ratio_med']=float(np.median(st['ratios'])) if st['ratios'] else 0; st['ratio_max']=float(np.max(st['ratios'])) if st['ratios'] else 0; st.pop('ratios')
    STATS[key]=st
    return st
for k in ('E1','E2','E3'):
    s=make(k,None); print(k, s)
json.dump(STATS, open(OUTDIR / 'mesh_stats.json', 'w', encoding='utf-8'), indent=1)

"""Reconcile every row of the analysis matrix against its component sheets.

`data/raw/mainfort-pfg-cpl.xlsx` sheet `pfg-cpl-mainfort` is the 55-assemblage
decorated-class matrix every transmission result rests on. It was built by
merging three sources, and nothing in the repository recorded which rows came
from where, so the provenance is measured here rather than described.

The answer, 2026-09-19: 35 rows are Mainfort (1996) alone, 19 sum Mainfort with
Lipo's compilation class by class, 1 is Lipo alone. In the 29-assemblage basin
set the split is 13 / 15 / 1. Nothing is left over.

The merge raises sample sizes; counts are summed, not chosen between. Seven of
Lipo's assemblages carry 1996-97 controlled surface collections and limited
excavation (Lipo 2001, Ch. 5), so part of what is summed is new material rather
than a recount of the same sherds.

The work is in the aliases. The three sheets key the same assemblage three
ways: `Kent_Place` in the merged sheet is `Kent` in Mainfort's and `13-N-4` in
Lipo's. Without them three rows look unreconcilable and six more look like
Mainfort-only rows.

Run: .venv/bin/python scripts/reconcile_ceramic_matrix.py
"""
import pandas as pd, numpy as np, re, itertools
from collections import Counter
def norm(s): return re.sub(r'[^a-z0-9]','',str(s).lower())
TY=['Parkin_Punctated','Barton/Kent/MPI','Painted','Fortune_Noded','Ranch_Incised','Walls_Engraved','Wallace_Incised','Rhodes_Incised','Vernon_Paul_Applique','Hull_Engraved']
F='data/raw/mainfort-pfg-cpl.xlsx'
ren={'ranch_inc':'Ranch_Incised','painted':'Painted','rhodes_inc':'Rhodes_Incised','walls_eng':'Walls_Engraved','vernon_pal':'Vernon_Paul_Applique','fortune_nd':'Fortune_Noded'}
# Name/grid aliases: the merged sheet, Mainfort's sheet and Lipo's sheet key the
# same assemblage three ways (Kent_Place / Kent / 13-N-4).
M_ALIAS={'kentplace':'Kent','bartonranch':'Barton','cramorplace':'Cramor'}
C_GRID={'kentplace':'13-N-4','lakecormorant':'13-P-8','castilelanding':'13-N-21','nickel':'13-N-15','rosemound':'12-N-3','holdenlake':'Holden_Lake'}
u=pd.read_excel(F,sheet_name='pfg-cpl-mainfort').dropna(subset=['Assemblages']); u['k']=u['Assemblages'].map(norm)
m=pd.read_excel(F,sheet_name='mainfort-data-collapsed').rename(columns=ren); m['n']=m['site_name'].astype(str).str.strip(); m['k']=m['n'].map(norm)
c=pd.read_excel(F,sheet_name='cpl-pfg'); c['id']=c['Assemblages'].astype(str).str.strip()
pfg=pd.read_excel(F,sheet_name='PFG'); pfg['Site Number']=pfg['Site Number'].astype(str).str.strip()
g2n={r['Site Number']:str(r['Site Name']).strip() for _,r in pfg.dropna(subset=['Site Name']).iterrows()}
def mainfort(k):
    name=M_ALIAS.get(k)
    sub=m[m['n']==name] if name else m[m['k']==k]
    return sub[[t for t in TY if t in m.columns]].sum().reindex(TY).fillna(0).values
def cpl(k):
    gid=C_GRID.get(k)
    ids=[gid] if gid else [i for i in c['id'] if norm(g2n.get(i,i))==k]
    sub=c[c['id'].isin(ids)]
    return sub[TY].apply(pd.to_numeric,errors='coerce').fillna(0).sum().reindex(TY).fillna(0).values
tal={}
for _,r in u.iterrows():
    k=r['k']; uv=r[TY].astype(float).values; src={'Mainfort':mainfort(k),'CPL':cpl(k)}
    hit='unreconciled'
    for n in (1,2):
        for combo in itertools.combinations(src,n):
            if np.allclose(sum(src[s] for s in combo),uv): hit='+'.join(combo); break
        if hit!='unreconciled': break
    tal[str(r['Assemblages'])]=hit
mem={ln.strip() for ln in open('data/processed/basin_members_curated.txt') if ln.strip()}
print("all 55:",dict(Counter(tal.values())))
print("basin 29:",dict(Counter(v for a,v in tal.items() if a in mem)))
print("unreconciled:",[a for a,v in tal.items() if v=='unreconciled'])

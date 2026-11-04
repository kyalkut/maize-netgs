import sys, numpy as np, csv
import os; HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); sys.path.insert(0,HERE)
from sbml import read_model
from netgs_step1 import genotype_fluxes
R=os.environ.get('NETGS_REPO','netGS/netGS/')
m=read_model(R+'model.xml'); S=m['S']; lb=m['lb']; ub=m['ub']; n=S.shape[1]; jb=548
bofs=np.loadtxt(R+'biomreaction.csv',delimiter=','); biom=np.loadtxt(R+'biomass_optN.csv',delimiter=',')[:,1]
carb,oxy,starch,suc=5,84,21,42   # 0-based of MATLAB 6,85,22,43
def row(i,j,ci,cj): r=np.zeros(n); r[i]=ci; r[j]=cj; return r
# reference: fixed ratios (equalities) in MATLAB step 1; relaxed ratio bounds for genotypes. Use relaxed bounds throughout here:
A=np.vstack([row(carb,oxy,1,-3.81445), -row(carb,oxy,1,-0.93815), row(starch,suc,1,-3.3694), -row(starch,suc,1,-0.7898)]); b=np.zeros(4)
N=int(sys.argv[1]) if len(sys.argv)>1 else 67
SOLVER=os.environ.get('NETGS_SOLVER','cvxpy')   # 'cvxpy' (Clarabel, like CVX), 'osqp' or 'highs' (code objective as LP)
res={}
for obj in ('paper','code'):
    r=genotype_fluxes(S,lb,ub,jb,bofs[:,:N],biom[:N],ref=14,extra_ub=(A,b),objective=obj,ref_bof=bofs[:,14],solver=SOLVER)
    res[obj]=r
    F=r['fluxes']; vr=r['vref'][r['nonzero']]; rel=F/vr[:,None]
    kb=int(np.where(r['nonzero']==jb)[0][0])
    print(f"[{obj}] nonzero reactions {len(r['nonzero'])} (MATLAB code: 336) | s_model {r['s_model']:.5g} | status {set(r['status'])}")
    print(f"   sum_k v_k/vref_k per genotype: median {np.nanmedian(np.nansum(rel,0)):.3g} (reference = {len(vr)})")
    print(f"   per-reaction relative deviation |v/vref-1|: median {np.nanmedian(np.abs(rel-1)):.3g}, 90th pct {np.nanpercentile(np.abs(rel-1),90):.3g}")
    print(f"   sign flips vs reference: {np.nanmean(rel<0)*100:.1f}% of reaction-genotype pairs")
    if obj=='code': print(f"   objective |sum v/vref - 1| per genotype: median {np.nanmedian(np.abs(np.nansum(rel,0)-1)):.6g}")
    print(f"   biomass flux vs target corr: {np.corrcoef(F[kb],r['target'][:N])[0,1]:.4f}")
sfx='' if SOLVER=='cvxpy' else '_'+SOLVER
np.save(os.path.join(HERE,f'arab_fluxes_paper{sfx}.npy'),res['paper']['fluxes']); np.save(os.path.join(HERE,f'arab_fluxes_code{sfx}.npy'),res['code']['fluxes'])

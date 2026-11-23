import sys, numpy as np
import os; HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); sys.path.insert(0,HERE)
import netgs_step1 as ng
from scipy.optimize import linprog as _lp
from sbml import read_model
R=os.environ.get('NETGS_REPO','netGS/netGS/')
m=read_model(R+'model.xml'); S=m['S']; n=S.shape[1]; jb=548
bofs=np.loadtxt(R+'biomreaction.csv',delimiter=','); biom=np.loadtxt(R+'biomass_optN.csv',delimiter=',')[:,1]
def row(i,j,ci,cj): r=np.zeros(n); r[i]=ci; r[j]=cj; return r
A=np.vstack([row(5,84,1,-3.81445),-row(5,84,1,-0.93815),row(21,42,1,-3.3694),-row(21,42,1,-0.7898)]); b=np.zeros(4)
P=np.load(os.path.join(HERE,'arab_fluxes_paper.npy')); C=np.load(os.path.join(HERE,'arab_fluxes_code.npy'))
# rerun 'code' objective with a different LP algorithm (interior point) -> uniqueness test
orig=ng.linprog
ng.linprog=lambda *a,**k: _lp(*a,**{**k,'method':'highs-ipm'})
r2=ng.genotype_fluxes(S,m['lb'],m['ub'],jb,bofs,biom,ref=14,extra_ub=(A,b),objective='code',ref_bof=bofs[:,14])
ng.linprog=orig
C2=r2['fluxes']
def per_reaction_corr(X,Y):
    out=[]
    for k in range(X.shape[0]):
        x,y=X[k],Y[k]; ok=np.isfinite(x)&np.isfinite(y)
        if np.std(x[ok])<1e-12 or np.std(y[ok])<1e-12: out.append(np.nan); continue
        out.append(np.corrcoef(x[ok],y[ok])[0,1])
    return np.array(out)
var=lambda X: np.mean(np.nanstd(X,1)>1e-9*np.nanmax(np.abs(X)))
print('reactions varying across genotypes: paper %.0f%%, code %.0f%%'%(100*var(P),100*var(C)))
pc=per_reaction_corr(P,C); cc=per_reaction_corr(C,C2)
print('per-reaction correlation across genotypes, paper vs code objective: median %.2f, share r>0.9: %.0f%%'%(np.nanmedian(pc),100*np.nanmean(pc>0.9)))
print('per-reaction correlation, code objective solved with two LP algorithms (simplex vs interior point): median %.2f, share r>0.9: %.0f%%'%(np.nanmedian(cc),100*np.nanmean(cc>0.9)))
print('max abs difference between the two code solutions (relative to flux scale): %.3g'%(np.nanmax(np.abs(C-C2))/np.nanmax(np.abs(C))))
np.save(os.path.join(HERE,'arab_fluxes_code_ipm.npy'),C2)

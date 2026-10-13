from pathlib import Path
import csv,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from figure_style_final import two_panel,style,save
ROOT=Path(__file__).resolve().parents[1]

def read_csv(name):
    with (ROOT/'data'/name).open(newline='') as f:return list(csv.DictReader(f))

def main():
    iso=read_csv('isolated_compacton_validation.csv')
    con=read_csv('absolute_spatial_convergence.csv')
    h=np.array([float(r['h']) for r in iso]); e3=np.array([float(r['relL2_c3']) for r in iso]); e1=np.array([float(r['relL2_c1']) for r in iso])
    nf=np.array([float(r['N_fine']) for r in con]); ef=np.array([float(r['registered_relL2']) for r in con]); el=np.array([float(r['low32_registered']) for r in con]); ec=np.array([float(r['E_core']) for r in con])
    pfield=np.log2(ef[0]/ef[1]);plow=np.log2(el[0]/el[1]);pcore=np.log2(ec[0]/ec[1])
    fig,ax1,ax2=two_panel()
    ax1.loglog(h,e3,'o-',label=r'$c=3$')
    ax1.loglog(h,e1,'s-',label=r'$c=1$')
    href=np.array([h.min(),h.max()]); cref=e1[-1]/h[-1]**2
    ax1.loglog(href,cref*href**2,'--',color='0.45',label=r'$h^2$ guide')
    ax1.invert_xaxis(); style(ax1,r'mesh size $h$',r'relative $L^2$ error','(a) Isolated compactons')
    ax1.legend(frameon=False,loc='upper right',handlelength=2.0)

    ax2.semilogy(nf,ef,'o-',label='registered field')
    ax2.semilogy(nf,el,'s-',label=r'low modes $|m|\leq32$')
    ax2.semilogy(nf,ec,'^-',label='fitted cores')
    ax2.set_xticks([1600,3200]); ax2.set_xticklabels(['1600','3200'])
    ax2.set_xlim(1350,3450); ax2.set_ylim(3.5e-4,1.5e-2)
    style(ax2,r'fine-grid size $N$',r'pairwise relative difference','(b) Absolute-flow collision')
    ax2.legend(frameon=False,loc='upper right',handlelength=1.8)
    save(fig,ROOT/'figures'/'fig5_absolute_spatial_convergence.pdf')
if __name__=='__main__':main()

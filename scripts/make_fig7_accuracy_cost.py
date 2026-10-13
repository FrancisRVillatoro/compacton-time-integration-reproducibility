from pathlib import Path
import csv,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from figure_style_final import one_panel,style,save
ROOT=Path(__file__).resolve().parents[1]

def main():
    with (ROOT/'data'/'table3_final.csv').open(newline='') as f: rows=list(csv.DictReader(f))
    methods=[('M2','o'),('RadauIIA2','s'),('TR-BDF2','^')]
    labels={'M2':'M2','RadauIIA2':'Radau IIA','TR-BDF2':'TR-BDF2'}
    fig,ax=one_panel()
    for method,marker in methods:
        rr=[r for r in rows if r['method']==method]
        rr=sorted(rr,key=lambda r:float(r['dt']),reverse=True)
        x=[float(r['wall_s']) for r in rr]; y=[1000*float(r['E_core']) for r in rr]
        ax.plot(x,y,marker=marker,label=labels[method])
    style(ax,'wall time [s]',r'$10^3 E_{core}$','Accuracy-cost comparison')
    ax.set_xlim(2.5,15.7);ax.set_ylim(1.9,12.0)
    ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1.02,1.0),borderaxespad=0.0,handlelength=2.0)
    save(fig,ROOT/'figures'/'fig7_accuracy_cost_absolute_reference.pdf')
if __name__=='__main__':main()

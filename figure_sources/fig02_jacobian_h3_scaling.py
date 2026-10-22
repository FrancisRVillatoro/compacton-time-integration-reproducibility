"""Figure 2, nonlinear Jacobian mesh scaling and stiff-mode amplification."""
from __future__ import annotations
import numpy as np
from plot_support import rows,column,single_panel,output


def main():
    allrows=rows('fig2_jacobian_h3_scaling.csv')
    sep=sorted((r for r in allrows if r['state']=='separated'),key=lambda r:float(r['h']))
    ov=sorted((r for r in allrows if r['state']=='overlap'),key=lambda r:float(r['h']))
    h=column(sep,'h')
    lam_s=column(sep,'abs_lambda_max');lam_o=column(ov,'abs_lambda_max')
    # The three finest meshes define both fitted exponents.
    slope_s=-np.polyfit(np.log(h[:3]),np.log(lam_s[:3]),1)[0]
    slope_o=-np.polyfit(np.log(h[:3]),np.log(lam_o[:3]),1)[0]
    fig,ax=single_panel()
    ax.loglog(h,lam_s,'s-',mfc='white',label=fr'separated, fit $h^{{-{slope_s:.2f}}}$')
    ax.loglog(h,lam_o,'D-',mfc='white',label=fr'overlap, fit $h^{{-{slope_o:.2f}}}$')
    href=np.array([h.min(),h.max()])
    ax.loglog(href,(lam_s[0]*.60)*(href/h[0])**(-3),'--',color='0.5')
    ax.text(.105,1.30e4,r'$\propto h^{-3}$',color='0.40',fontsize=9)
    ax.set_xlabel(r'mesh size $h$');ax.set_ylabel(r'$|\lambda_{\max}|$')
    ax.set_xticks([.025,.05,.1,.2,.4]);ax.set_xticklabels(['0.025','0.05','0.1','0.2','0.4'])
    ax.set_xlim(.022,.46);ax.set_ylim(350,1.4e7)
    ax.legend(frameon=False,loc='upper right',fontsize=8.1,handlelength=1.3)
    output(fig,'fig2_jacobian_h3_scaling_left.pdf')

    fig,ax=single_panel()
    series=[('M3_amp','^','C0','fourth-order split-step'),
            ('M1_amp','o','C1','implicit midpoint'),
            ('GA_rho05_amp','o','C2',r'gen.-$\alpha$, $\rho_\infty=0.5$'),
            ('GA_rho0_amp','p','C3',r'gen.-$\alpha$, $\rho_\infty=0$'),
            ('M2_amp','s','C4',r'$R_{1,2}$ Padé')]
    for key,marker,color,label in series:
        ax.loglog(h,column(sep,key),marker+'-',color=color,label=label,markersize=5)
    # Reference decay slopes, not fitted confidence intervals.
    ax.loglog(href,.07*(href/.4)**3,'--',color='0.5',lw=1.0)
    ax.loglog(href,.10*(href/.4)**1.5,'--',color='0.5',lw=1.0)
    ax.text(.11,3.e-5,r'$\propto h^{3}$',fontsize=9,color='0.4')
    ax.text(.095,7.e-3,r'$\propto h^{3/2}$',fontsize=9,color='0.4')
    ax.set_xlabel(r'mesh size $h$');ax.set_ylabel('stiff-mode amplification')
    ax.set_xticks([.025,.05,.1,.2,.4]);ax.set_xticklabels(['0.025','0.05','0.1','0.2','0.4'])
    ax.set_xlim(.022,.46);ax.set_ylim(7.e-6,2)
    ax.legend(frameon=False,loc='lower right',fontsize=8,handlelength=1.5,labelspacing=.2)
    output(fig,'fig2_jacobian_h3_scaling_right.pdf')


if __name__=='__main__': main()

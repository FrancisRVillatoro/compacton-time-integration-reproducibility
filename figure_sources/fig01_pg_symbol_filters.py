"""Figure 1, exact PG dispersion and method amplification (two PDF panels)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from plot_support import rows, column, single_panel, output

# Figure 1 and the Jacobian experiment share one canonical generalized-alpha
# amplification matrix.  No nonlinear generalized-alpha stepper is used here.
JACOBIAN_SOURCE = Path(__file__).resolve().parents[1] / "reproducibility" / "jacobian"
sys.path.insert(0, str(JACOBIAN_SOURCE))
from jacobian_core import ga_rho  # noqa: E402


def ga_amplification(z: np.ndarray, rho_inf: float) -> np.ndarray:
    return np.array([ga_rho(complex(value), rho_inf) for value in z], dtype=float)


def main():
    data=rows('fig1_pg_symbol_filters.csv')
    x=column(data,'theta_over_pi')
    zim=column(data,'z_q5625_imag')/1.e4
    cont=column(data,'z_cont_imag')/1.e4
    q0=column(data,'z_q0_imag')/1.e4
    fig,ax=single_panel()
    ax.plot(x,zim,color='C0',label=r'exact PG symbol, $q=5.625$')
    ax.plot(x,cont,'--',color='C1',label=r'continuum proxy, $q=5.625$')
    ax.plot(x,q0,':',color='0.40',label=r'exact PG symbol, $q=0$')
    i=int(np.argmax(zim))
    ax.plot(x[i],zim[i],marker='o',markerfacecolor='none',markeredgecolor='C2',markersize=4.5,zorder=5)
    ax.annotate(r'$\theta^*\simeq0.80\pi$',xy=(x[i],zim[i]),xytext=(.62,7.1),
                arrowprops=dict(arrowstyle='-',lw=.8,color='0.20'),fontsize=9,ha='left')
    ax.set_xlim(0,1);ax.set_ylim(-.45,14.5)
    ax.legend(frameon=False,loc='upper left',handlelength=2.0,fontsize=8.4)
    ax.set_xlabel(r'$\theta/\pi$');ax.set_ylabel(r'$\operatorname{Im}z(\theta)/10^4$')
    output(fig,'fig1_pg_symbol_filters_left.pdf')

    # The CSV contains the exact PG scalar amplification of midpoint, Padé,
    # Radau (identical to Padé), TR-BDF2, and the Padé continuum comparison.
    z=column(data,'z_q5625_real') + 1j*column(data,'z_q5625_imag')
    fig,ax=single_panel()
    ax.semilogy(x,column(data,'amp_midpoint_exact'),color='C0',label='midpoint, split-step')
    ax.semilogy(x,ga_amplification(z,.5),color='C1',label=r'gen.-$\alpha$, $\rho_\infty=0.5$')
    ax.semilogy(x,ga_amplification(z,0.),color='C2',label=r'gen.-$\alpha$, $\rho_\infty=0$')
    ax.semilogy(x,column(data,'amp_TRBDF2_exact'),color='C3',label='TR-BDF2')
    ax.semilogy(x,column(data,'amp_M2_Radau_exact'),color='C4',label=r'$R_{1,2}$ Padé, Radau IIA')
    ax.semilogy(x,column(data,'amp_M2_continuum'),'--',color='0.50',label=r'$R_{1,2}$ Padé, continuum proxy')
    ax.set_xlabel(r'$\theta/\pi$');ax.set_ylabel(r'amplification $|R(z(\theta))|$')
    ax.set_xlim(0,1);ax.set_ylim(1.e-5,1.8)
    ax.legend(frameon=False,loc='center right',fontsize=8.0,handlelength=1.65,labelspacing=.2)
    output(fig,'fig1_pg_symbol_filters_right.pdf')


if __name__=='__main__':main()

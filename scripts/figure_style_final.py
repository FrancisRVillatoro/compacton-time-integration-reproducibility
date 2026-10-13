from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timezone
MM=1/25.4
AX_W=72*MM
AX_H=54*MM
DOUBLE_W=190*MM
ONEHALF_W=140*MM
BOTTOM=13*MM
TOP=5*MM
LEFT=14*MM
GAP=23*MM
H=BOTTOM+AX_H+TOP
FONT=10
mpl.rcParams.update({
 'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
 'font.size':FONT,'axes.labelsize':FONT,'axes.titlesize':FONT,'xtick.labelsize':FONT,
 'ytick.labelsize':FONT,'legend.fontsize':FONT,'pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':0.9,'lines.linewidth':1.45,'lines.markersize':6.5,
 'figure.facecolor':'white','axes.facecolor':'white','savefig.bbox':'standard',
 'mathtext.default':'regular'
})
def style(ax,xlabel=None,ylabel=None,title=None):
    ax.minorticks_on(); ax.tick_params(which='major',direction='out',length=4,width=.9)
    ax.tick_params(which='minor',direction='out',length=2.5,width=.7)
    for s in ax.spines.values(): s.set_linewidth(.9)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if title: ax.set_title(title,pad=6)
    return ax
def two_panel():
    fig=plt.figure(figsize=(DOUBLE_W,H))
    ax1=fig.add_axes([LEFT/DOUBLE_W,BOTTOM/H,AX_W/DOUBLE_W,AX_H/H])
    ax2left=LEFT+AX_W+GAP
    ax2=fig.add_axes([ax2left/DOUBLE_W,BOTTOM/H,AX_W/DOUBLE_W,AX_H/H])
    return fig,style(ax1),style(ax2)
def one_panel():
    fig=plt.figure(figsize=(ONEHALF_W,H))
    ax=fig.add_axes([LEFT/ONEHALF_W,BOTTOM/H,AX_W/ONEHALF_W,AX_H/H])
    return fig,style(ax)
def save(fig,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    meta={'CreationDate':datetime(2026,10,10,12,0,0,tzinfo=timezone.utc),'ModDate':datetime(2026,10,10,12,0,0,tzinfo=timezone.utc)}
    fig.savefig(path,facecolor='white',metadata=meta);plt.close(fig)

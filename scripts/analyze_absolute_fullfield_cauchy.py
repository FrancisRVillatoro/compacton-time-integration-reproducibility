#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.signal import resample
from scipy.optimize import minimize_scalar


def load_npz(path):
    z=np.load(path,allow_pickle=True)
    meta=json.loads(str(z['metadata'].item()))
    return dict(path=str(path), x=z['x'], times=z['times'], snapshots=z['snapshots'], meta=meta)


def periodic_shift(u, shift, L):
    N=u.size
    k=2*np.pi*np.fft.fftfreq(N, d=L/N)
    U=np.fft.fft(u)
    return np.fft.ifft(U*np.exp(-1j*k*shift)).real


def rel(a,b):
    den=np.linalg.norm(b)
    return float(np.linalg.norm(a-b)/den)


def best_shift(a,b,L):
    # Shift a to match b. Coarse estimate by circular correlation, then scalar refinement.
    N=b.size
    Fa=np.fft.fft(a); Fb=np.fft.fft(b)
    corr=np.fft.ifft(Fb*np.conj(Fa)).real
    j=int(np.argmax(corr))
    if j>N//2: j-=N
    dx=L/N
    s0=j*dx
    # Search a narrow local interval; objective is smooth under Fourier shift.
    f=lambda s: rel(periodic_shift(a,s,L),b)
    opt=minimize_scalar(f, bounds=(s0-2*dx,s0+2*dx), method='bounded', options={'xatol':1e-12})
    return float(opt.x), float(opt.fun)


def lowpass(u,mcut):
    N=u.size
    U=np.fft.fft(u)
    m=np.fft.fftfreq(N)*N
    U[np.abs(m)>mcut]=0
    return np.fft.ifft(U).real


def neg_frac(u):
    return float(np.linalg.norm(np.minimum(u,0.0))/np.linalg.norm(u))


def pos_frac(u):
    return float(np.linalg.norm(np.maximum(u,0.0))/np.linalg.norm(u))


def compare(coarse,fine,selected_times=(0,2,3.7,4,5,6,7,7.5,7.75,8,8.5,9,9.5,10)):
    mc,mf=coarse['meta'],fine['meta']
    for key in ['L','T','c0','mu','cfl','eps']:
        if abs(float(mc[key])-float(mf[key]))>1e-13:
            raise RuntimeError(f'metadata mismatch {key}: {mc[key]} vs {mf[key]}')
    L=float(mf['L'])
    tc=coarse['times']; tf=fine['times']
    if tc.shape!=tf.shape or np.max(np.abs(tc-tf))>1e-12:
        raise RuntimeError('time grids differ')
    Nf=fine['snapshots'].shape[1]
    hist=[]
    for i,t in enumerate(tf):
        a=resample(coarse['snapshots'][i],Nf)
        b=fine['snapshots'][i]
        raw=rel(a,b)
        sh,reg=best_shift(a,b,L)
        ar=periodic_shift(a,sh,L)
        row={'time':float(t),'raw_relL2':raw,'global_shift':sh,'registered_relL2':reg,
             'fine_negative_L2_fraction':neg_frac(b),'coarse_negative_L2_fraction':neg_frac(a),
             'fine_min_u':float(np.min(b)),'coarse_min_u':float(np.min(a))}
        for mcut in (16,32,64,128):
            al=lowpass(a,mcut); bl=lowpass(b,mcut)
            sl,er=best_shift(al,bl,L)
            row[f'low_m{mcut}_shift']=sl
            row[f'low_m{mcut}_registered']=er
        # Difference decomposition after global registration.
        row['positive_registered_rel_to_fullfine']=float(np.linalg.norm(np.maximum(ar,0)-np.maximum(b,0))/np.linalg.norm(b))
        row['negative_registered_rel_to_fullfine']=float(np.linalg.norm(np.minimum(ar,0)-np.minimum(b,0))/np.linalg.norm(b))
        hist.append(row)
    # nearest selected times
    selected=[]
    for ts in selected_times:
        j=int(np.argmin(np.abs(tf-ts)))
        selected.append(hist[j])
    return hist,selected


def write_csv(rows,path):
    import csv
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def summarize(label,coarse,fine,hist,selected,core_pair=None):
    last=hist[-1]
    lines=[]
    lines.append(f'{label}: N={coarse["meta"]["N"]} -> {fine["meta"]["N"]}, absolute |u|u, T={fine["meta"]["T"]}')
    lines.append(f'final raw_relL2={last["raw_relL2"]:.12e}')
    lines.append(f'final registered_relL2={last["registered_relL2"]:.12e} global_shift={last["global_shift"]:+.12e}')
    for m in (16,32,64,128):
        lines.append(f'final low_m{m}_registered={last[f"low_m{m}_registered"]:.12e} shift={last[f"low_m{m}_shift"]:+.12e}')
    lines.append(f'final positive_difference/fullfine={last["positive_registered_rel_to_fullfine"]:.12e}')
    lines.append(f'final negative_difference/fullfine={last["negative_registered_rel_to_fullfine"]:.12e}')
    lines.append(f'final negative_L2_fraction coarse={last["coarse_negative_L2_fraction"]:.12e} fine={last["fine_negative_L2_fraction"]:.12e}')
    lines.append(f'final min_u coarse={last["coarse_min_u"]:.12e} fine={last["fine_min_u"]:.12e}')
    if core_pair is not None:
        lines.append(f'core_pair_difference={core_pair:.12e}')
        lines.append(f'full_registered/core_pair={last["registered_relL2"]/core_pair:.6f}')
    lines.append('SELECTED TIMES')
    for r in selected:
        lines.append(f't={r["time"]:.3f} raw={r["raw_relL2"]:.6e} reg={r["registered_relL2"]:.6e} low32={r["low_m32_registered"]:.6e} negFine={r["fine_negative_L2_fraction"]:.6e} minFine={r["fine_min_u"]:.6e}')
    return '\n'.join(lines)+'\n'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('coarse')
    ap.add_argument('fine')
    ap.add_argument('--out-prefix',default='absolute_fullfield')
    ap.add_argument('--core-pair',type=float,default=None)
    args=ap.parse_args()
    c=load_npz(args.coarse); f=load_npz(args.fine)
    if int(c['meta'].get('ext',-1))!=1 or int(f['meta'].get('ext',-1))!=1:
        raise SystemExit('both inputs must be ext=1 absolute_K22_localization')
    hist,sel=compare(c,f)
    prefix=Path(args.out_prefix)
    write_csv(hist,str(prefix)+'_history.csv')
    report=summarize(prefix.name,c,f,hist,sel,args.core_pair)
    Path(str(prefix)+'_report.txt').write_text(report)
    print(report)

if __name__=='__main__':
    main()

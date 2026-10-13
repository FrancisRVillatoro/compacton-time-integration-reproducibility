#!/usr/bin/env python3
"""Audit the I2 balance along stored PG states; does NOT integrate the PDE.

Required NPZ arrays: times, U (or snapshots), x, metadata (JSON scalar).
Metadata: L, p, c0, mu, N. No domain/model parameters are guessed.
The trapezoidal correction is a sampled balance diagnostic, not a certified
pure time-discretization error. Coarsening and Simpson comparisons are
reported so that undersampled defects cannot silently pass a threshold.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib, json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson


def trapz(y, t):
    return float(np.sum(0.5*(y[:-1]+y[1:])*np.diff(t)))


def cumulative_trapz(y, t):
    return np.r_[0.0, np.cumsum(0.5*(y[:-1]+y[1:])*np.diff(t))]


def load_model_module(source_dir):
    root=Path(source_dir).resolve()
    expected=root/'pg_lstable_methods.py'
    if not expected.is_file():
        raise FileNotFoundError(expected)
    sys.path.insert(0,str(root))
    pm=importlib.import_module('pg_lstable_methods')
    if Path(pm.__file__).resolve()!=expected:
        raise RuntimeError('An unexpected pg_lstable_methods module was imported')
    return pm


def diagnostic_series(model,U):
    h,p,c0,mu=model.h,model.p,model.c0,model.mu
    vals=[]
    for u in U:
        g=model.g(u)
        i2=float(h*np.sum(np.abs(u)**(p+1)/(p+1)))
        adv=float(c0*h*np.dot(g,model.Alu.solve(model.B@u)))
        visc=float(-mu*h**3*np.dot(g,model.Alu.solve(model.D@u))) if mu else 0.0
        nonlinear=float(-h*np.dot(g,model.Alu.solve(model.BC@g)))
        direct=float(-h*np.dot(g,model.F(u)))
        vals.append([i2,adv,visc,nonlinear,direct])
    return np.asarray(vals)


def audit(path,source_dir=None,out_prefix=None):
    path=Path(path)
    with np.load(path,allow_pickle=False) as z:
        if 'times' not in z or not ('U' in z or 'snapshots' in z):
            raise ValueError('This archive does not contain a trajectory: initial/final states alone are insufficient.')
        t=np.asarray(z['times'],float)
        U=np.asarray(z['U'] if 'U' in z else z['snapshots'],float)
        x=np.asarray(z['x'],float)
        if 'metadata' not in z:
            raise ValueError('Missing explicit model metadata')
        meta=json.loads(str(z['metadata'].item()))
    for key in ('L','p','c0','mu','N'):
        if key not in meta: raise ValueError('Missing parameter: '+key)
    N=int(meta['N']); L=float(meta['L']); h=L/N
    if U.ndim!=2 or U.shape!=(len(t),N) or x.shape!=(N,):
        raise ValueError('Unexpected state/grid dimensions')
    if len(t)<3 or np.any(np.diff(t)<=0) or not np.all(np.isfinite(U)):
        raise ValueError('Need at least three ordered finite states')
    if not np.allclose(np.diff(x),h,rtol=1e-11,atol=1e-12):
        raise ValueError('Grid spacing is inconsistent with L/N')
    if abs(t[0])>1e-13 or ('T' in meta and abs(t[-1]-meta['T'])>1e-10):
        raise ValueError('Trajectory must include initial and final times')
    source_dir=Path(source_dir or Path(__file__).parent)
    pm=load_model_module(source_dir)
    model=pm.AbsKPetrovGalerkin(N,L,float(meta['p']),float(meta['c0']),float(meta['mu']))
    vals=diagnostic_series(model,U)
    I20=vals[0,0]
    if I20<=0: raise ValueError('I2(0) must be positive')
    adv=vals[:,1]/I20; visc=vals[:,2]/I20; semi=adv+visc
    raw=(vals[:,0]-I20)/I20
    integ=cumulative_trapz(semi,t)
    corrected=raw-integ
    # Independent FFT evaluation for this circulant PG stencil.
    theta=2*np.pi*np.fft.fftfreq(N)
    a=(66+52*np.cos(theta)+2*np.cos(2*theta))/120
    b=1j*(10*np.sin(theta)+np.sin(2*theta))/(12*h)
    fft_adv=np.asarray([model.c0*h*np.dot(model.g(u),np.fft.ifft((b/a)*np.fft.fft(u)).real)/I20 for u in U])
    sampling=[]
    for stride in (1,2,4,5,10,20):
        ind=np.unique(np.r_[np.arange(0,len(t),stride),len(t)-1])
        val=trapz(semi[ind],t[ind])
        sampling.append(dict(stride=stride,samples=len(ind),integral_signed_relative=val,
                             abs_integral_below_1e_8=bool(abs(val)<1e-8)))
    main_int=float(integ[-1]); simp=float(simpson(semi,x=t))
    report=dict(
        trajectory_file=path.name,trajectory_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(pm.__file__).read_bytes()).hexdigest(),
        model_metadata=meta,states=len(t),I2_initial=I20,I2_final=float(vals[-1,0]),
        raw_drift_signed_relative=float(raw[-1]),
        advection_defect_trapezoid_signed_relative=trapz(adv,t),
        viscosity_defect_trapezoid_signed_relative=trapz(visc,t),
        semidiscrete_defect_trapezoid_signed_relative=main_int,
        corrected_drift_trapezoid_signed_relative=float(corrected[-1]),
        corrected_drift_trapezoid_absolute_relative=float(abs(corrected[-1])),
        absolute_defect_to_absolute_drift_ratio=abs(main_int)/abs(raw[-1]) if raw[-1] else None,
        max_absolute_cumulative_defect_relative=float(np.max(abs(integ))),
        integral_absolute_defect_relative=trapz(abs(semi),t),
        simpson_signed_relative=simp,simpson_minus_trapezoid=simp-main_int,
        sampling_checks=sampling,
        max_fft_lu_difference_relative=float(np.max(abs(fft_adv-adv))),
        max_nonlinear_cancellation_relative=float(np.max(abs(vals[:,3]))/I20),
        max_chain_rule_identity_residual_relative=float(np.max(abs(vals[:,4]-vals[:,1]-vals[:,2]))/I20),
        below_threshold_on_full_saved_grid=bool(abs(main_int)<1e-8),
        temporal_quadrature_certified=False,
        warning='A sampled trapezoidal balance is not a certified temporal-error separation. Inspect sampling_checks; do not treat cancellation after subsampling as an upper bound.'
    )
    prefix=Path(out_prefix) if out_prefix else path.with_name(path.stem+'_I2_defect')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    csv_path=prefix.with_suffix('.csv')
    with csv_path.open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['time','I2','delta_advection','delta_viscosity','nonlinear_cancellation','I2dot_direct','delta_advection_relative','delta_total_relative','I2_raw_drift_relative','cumulative_defect_relative','corrected_drift_relative'])
        for k in range(len(t)):
            w.writerow([t[k],*vals[k],adv[k],semi[k],raw[k],integ[k],corrected[k]])
    prefix.with_suffix('.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    text=[f"I2 DEFECT AUDIT: {path.name}",f"N={N} L={L:g} p={meta['p']} c0={meta['c0']} mu={meta['mu']} states={len(t)}",
          f"I2(0)                         {I20:.15e}",
          f"raw drift / I2(0)             {raw[-1]:+.15e}",
          f"signed integral delta/I2(0)   {main_int:+.15e}",
          f"corrected sampled drift       {corrected[-1]:+.15e}",
          f"Simpson                       {simp:+.15e}",
          'SAMPLING (not an error bound):']
    for r in sampling:
        text.append(f"  stride={r['stride']:2d} samples={r['samples']:4d} Q={r['integral_signed_relative']:+.15e} below_1e-8={r['abs_integral_below_1e_8']}")
    text.extend([report['warning'],'I2_DEFECT_AUDIT_COMPUTED'])
    prefix.with_suffix('.txt').write_text('\n'.join(text)+'\n')
    return report


def selftest(source_dir):
    pm=load_model_module(source_dir)
    rng=np.random.default_rng(1847)
    model=pm.AbsKPetrovGalerkin(64,100,2,1,0)
    u=rng.normal(size=64)
    v=diagnostic_series(model,np.array([u,np.ones(64)]))
    assert abs(v[0,4]-v[0,1])<2e-12*max(1,abs(v[0,1]))
    assert np.max(abs(v[1,1:]))<1e-13
    assert np.allclose(model.g(np.array([-2.,2.])),[-4,4])
    assert trapz(np.array([1.,3.,5.]),np.array([0.,1.,2.]))==6.
    print('I2_DEFECT_SELFTEST_OK')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('trajectory',nargs='?',type=Path)
    ap.add_argument('--source-dir',type=Path,default=Path(__file__).parent)
    ap.add_argument('--out-prefix',type=Path)
    ap.add_argument('--selftest',action='store_true')
    args=ap.parse_args()
    if args.selftest:
        selftest(args.source_dir)
    if args.trajectory is not None:
        result=audit(args.trajectory,args.source_dir,args.out_prefix)
        print(json.dumps(result,indent=2,allow_nan=False))
    elif not args.selftest:
        ap.error('provide a trajectory or --selftest')

if __name__=='__main__': main()

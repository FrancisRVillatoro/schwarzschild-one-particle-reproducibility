#!/usr/bin/env python3
"""Paper-III reproducibility runner.

Runs ALL scripts in an isolated temporary checkout, so frozen baselines and
SHA256SUMS are not rewritten during verification. Tests are not a proof of
any infinite-dimensional analytical theorem.
"""
from pathlib import Path
import argparse, csv, json, math, shutil, subprocess, sys, tempfile, time

PAPER=Path(__file__).resolve().parent
QUICK=[
    ('orientation', 'work1/check_time_orientation.py'),
    ('symbols', 'work2/checks/check_symbol_algebra.py'),
    ('spin_symbolic', 'work3/checks/check_work3_symbolics.py'),
    ('CAR_global_local', 'work8/checks/check_priority4.py'),
    ('massive_transition', 'work7/checks_work7.py'),
]
FULL=QUICK+[
    ('invariant_harmonics', 'work3/checks/check_invariant_harmonics.py'),
    ('continuum_packet', 'work3/checks/check_continuum_stress.py'),
    ('continuum_refinements', 'work3/checks/check_continuum_refinements.py'),
    ('Ward_CAR_continuum', 'work4/checks/check_work4.py'),
    ('transition_F_refinements', 'work5/checks/F_refinements.py'),
    ('transition_C_F_quadrature', 'work5/checks/convergence_CF.py'),
    ('transition_C_F_tail', 'work5/checks/CF_tail_refinement.py'),
    ('finite_window_coefficient', 'work6/checks/window_CF.py'),
    ('HS_toy_multiplicity', 'work6/checks/check_hs_counting.py'),
]

def check_frozen_baselines(p):
    """Check nonmonotone universal F, updated C_F and window coefficient."""
    results={}
    path=p/'work5/checks/F_refinements.csv'
    with path.open() as f: rows={float(r['Lambda']):r for r in csv.DictReader(f)}
    results['nonmonotonicity']=float(rows[4.4]['F_refined'])<float(rows[4.6]['F_refined'])<float(rows[5.0]['F_refined'])
    path=p/'work5/checks/CF_cutoff_convergence.csv'
    with path.open() as f: cf=list(csv.DictReader(f))
    results['CF_universal_model']=abs(float(cf[-1]['C_estimate'])-0.0788849362)<5e-8
    path=p/'work6/checks/window_integrals.csv'
    with path.open() as f: win=list(csv.DictReader(f))
    row=win[-1]
    results['window_positive']=float(row['I_window'])>0
    results['window_factor']=abs(float(row['8_over_3pi_times_I'])-float(row['I_window'])*8/(3*math.pi))<5e-12
    if not all(results.values()):raise AssertionError(f'Frozen reference outputs failed self-consistency: {results}')
    return results

def quick_f_check(p):
    """Fresh ODE computation at the non-monotone window, not only CSV read."""
    sys.path.insert(0,str(p/'work5/checks'))
    from universal_precision import transition
    vals={}
    for z in (4.4,4.6,5.0):
        r=transition(z,rtol=2e-11,order=2)
        vals[str(z)]={'F':r['F'],'norm_error':r['norm_error'],'nfev':r['nfev']}
        if abs(r['norm_error'])>5e-9:raise AssertionError(f'Unitary norm defect {z}: {r}')
    if not vals['4.4']['F']<vals['4.6']['F']<vals['5.0']['F']:
        raise AssertionError(f'Nonmonotonicity failed: {vals}')
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=('quick','full'),default='quick')
    ap.add_argument('--report',help='Optional explicit path for JSON report outside the repository')
    a=ap.parse_args()
    started=time.perf_counter()
    summary=dict(mode=a.mode,python=sys.version.split()[0],base='Paper III v1.2.0 candidate',checks={},steps=[],comments='Scripts are run in a disposable copy, not against the released files.')
    with tempfile.TemporaryDirectory(prefix='fermions_paperIII_') as t:
        p=Path(t)/'paper_III'
        shutil.copytree(PAPER,p)
        summary['checks']=check_frozen_baselines(p)
        todo=FULL if a.mode=='full' else QUICK
        for i,(label,name) in enumerate(todo,1):
            st=time.perf_counter()
            print(f'[{i}/{len(todo)}] {label}',flush=True)
            path=p/name
            proc=subprocess.run([sys.executable,str(path)],cwd=path.parent,capture_output=True,text=True)
            elapsed=time.perf_counter()-st
            details={'name':label,'file':name,'returncode':proc.returncode,'elapsed_s':round(elapsed,3),'stdout_tail':proc.stdout[-1800:],'stderr_tail':proc.stderr[-800:]}
            summary['steps'].append(details)
            if proc.returncode:
                print(proc.stdout[-3000:],proc.stderr[-3000:],flush=True)
                raise RuntimeError(f'CHECK FAILED {label}: {proc.returncode}')
            print(f'   PASS ({elapsed:.1f}s)',flush=True)
        st=time.perf_counter()
        summary['fresh_transition_window']=quick_f_check(p)
        print('   Fresh transition window: PASS',round(time.perf_counter()-st,2),'s',flush=True)
    summary['elapsed_s']=round(time.perf_counter()-started,3)
    summary['passed']=True
    if a.report:
        r=Path(a.report).expanduser().resolve();r.parent.mkdir(parents=True,exist_ok=True)
        r.write_text(json.dumps(summary,indent=2)+'\n')
        print('Report:',r,flush=True)
    print(f'Paper III {a.mode.upper()} checks: PASS ({summary["elapsed_s"]:.1f} seconds)',flush=True)

if __name__=='__main__':main()

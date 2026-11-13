#!/usr/bin/env python3
"""Paper-III reproducibility runner.

Runs ALL scripts in an isolated temporary checkout, so frozen baselines and
SHA256SUMS are not rewritten during verification. Tests are not a proof of
any infinite-dimensional analytical theorem.
"""
from pathlib import Path
import argparse, csv, json, math, os, shutil, subprocess, sys, tempfile, time

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

GENERATED_TABLES=[
    'work3/checks/invariant_harmonics_results.csv',
    'work3/checks/continuum_packet_results.csv',
    'work3/checks/continuum_refinements.csv',
    'work4/checks/box_limit.csv',
    'work4/checks/continuum_charges.csv',
    'work4/checks/terminal_bounds.csv',
    'work4/checks/ward_identity_refinement.csv',
    'work5/checks/F_refinements.csv',
    'work5/checks/CF_quadrature_convergence.csv',
    'work5/checks/CF_cutoff_convergence.csv',
    'work6/checks/window_integrals.csv',
    'work6/checks/toy_HS_counting.csv',
]

TABLE_RTOL=5e-9
TABLE_ATOL=5e-11
TERMINAL_FIELDS=('mu','r','scaled_density','density_limit','scaled_T12','T12_limit',
                 'additive_error','explicit_error_bound','bound_ratio')
TERMINAL_DERIVED={'additive_error','bound_ratio'}
TERMINAL_RS=2.3
TERMINAL_OPERAND_ULPS=8
TERMINAL_RATIO_ULPS=2
TERMINAL_BOUND_SLACK=1e-7


def terminal_diagnostics(row, location):
    """Check the two derived columns using the Work 4 formulae."""
    try:
        values={key:float(value) for key,value in row.items()}
    except ValueError as exc:
        raise AssertionError(f'Nonnumeric terminal bound data: {location}') from exc
    if not all(math.isfinite(value) for value in values.values()):
        raise AssertionError(f'Nonfinite terminal bound data: {location}')
    r=values['r']
    d,b,q=(values[key] for key in ('additive_error','explicit_error_bound','bound_ratio'))
    if not (0<r<1 and d>=0 and b>0 and q>=0):
        raise AssertionError(f'Invalid terminal bound data: {location}')
    if not d<b*(1+TERMINAL_BOUND_SLACK):
        raise AssertionError(f'Terminal error exceeds bound: {location}: {d} vs {b}')

    # With R=Rs*r, the raw leading density is density_limit / scale:
    # q_geometry*R = sqrt(Rs)*scale, where scale=R**2.5*sqrt(1-r).
    # Reversing the scale and using this alternate expression for q_geometry
    # involves several rounded operations (up to 3 operand ULPs in the frozen
    # rows). Use 8 operand ULPs to cover those operations, before subtraction.
    scale=(TERMINAL_RS*r)**2.5*math.sqrt(1-r)
    if scale<=0:
        raise AssertionError(f'Invalid terminal density scale: {location}')
    density=values['scaled_density']/scale
    leading=values['density_limit']/scale
    reconstructed=abs(density-leading)
    if not math.isfinite(reconstructed):
        raise AssertionError(f'Nonfinite reconstructed terminal residual: {location}')
    operand_ulp=max(math.ulp(density),math.ulp(leading))
    residual_budget=TERMINAL_OPERAND_ULPS*operand_ulp
    if abs(d-reconstructed)>residual_budget:
        raise AssertionError(f'Inconsistent terminal additive_error: {location}: '
                             f'{d} vs {reconstructed} (roundoff budget {residual_budget})')

    reconstructed_ratio=d/b
    ratio_budget=TERMINAL_RATIO_ULPS*math.ulp(reconstructed_ratio)
    if abs(q-reconstructed_ratio)>ratio_budget:
        raise AssertionError(f'Inconsistent terminal bound_ratio: {location}: '
                             f'{q} vs {reconstructed_ratio} (roundoff budget {ratio_budget})')
    return {'additive_error':d,'explicit_error_bound':b,'bound_ratio':q,
            'reconstructed_additive_error':reconstructed,'operand_ulp':operand_ulp,
            'additive_roundoff_budget':residual_budget,'ratio_roundoff_budget':ratio_budget}


def compare_terminal_diagnostics(actual, reference, row_index):
    location=f'terminal_bounds.csv:{row_index}'
    u=terminal_diagnostics(actual,location+':actual')
    v=terminal_diagnostics(reference,location+':reference')
    du,dv=u['additive_error'],v['additive_error']
    bu,bv=u['explicit_error_bound'],v['explicit_error_bound']
    qu,qv=u['bound_ratio'],v['bound_ratio']
    d_tol=max(TABLE_ATOL,TABLE_RTOL*max(du,dv),
              u['additive_roundoff_budget']+v['additive_roundoff_budget'])
    # The bound must pass the ordinary comparison as well. Propagate its
    # observed delta and the permitted residual delta through the quotient.
    q_tol=max(TABLE_ATOL,TABLE_RTOL*max(qu,qv),
              (d_tol+max(qu,qv)*abs(bu-bv))/min(bu,bv)
              +u['ratio_roundoff_budget']+v['ratio_roundoff_budget'])
    for key,x,y,tolerance in (('additive_error',du,dv,d_tol),('bound_ratio',qu,qv,q_tol)):
        if abs(x-y)>tolerance:
            raise AssertionError(f'Numerical mismatch {location}:{key}: {x} != {y} '
                                 f'(allowed absolute difference {tolerance})')
    return {'row':row_index,'mu':float(actual['mu']),'r':float(actual['r']),
            'actual':u,'reference':v,'additive_error_abs_difference':abs(du-dv),
            'additive_error_tolerance':d_tol,'bound_ratio_abs_difference':abs(qu-qv),
            'bound_ratio_tolerance':q_tol}


def compare_generated_table(actual, reference):
    """Compare fresh numerical results, retaining the frozen files unchanged.

    The default tolerances match the combined scalar CSV validator. Small F and
    tail values use a tighter absolute floor, so a zeroed tail cannot pass.
    Wall-clock times and adaptive evaluation counts are performance metadata.
    Norm defects are checked against the same 5e-9 bound as quick_f_check,
    rather than requiring reproduction of round-off noise bit for bit. The two
    derived terminal-bound columns account for roundoff before subtraction.
    """
    def read(path):
        with path.open(newline='') as stream:
            reader=csv.DictReader(stream);rows=list(reader);header=reader.fieldnames
        if not header or len(header)!=len(set(header)) or not rows:
            raise AssertionError(f'Empty or malformed table: {path}')
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise AssertionError(f'Malformed row: {path}')
        return header,rows
    ah,a=read(actual);bh,b=read(reference)
    if ah!=bh or len(a)!=len(b):raise AssertionError(f'Table shape differs: {actual.name}')
    terminal=actual.name==reference.name=='terminal_bounds.csv'
    terminal_rows=[]
    if terminal and tuple(ah)!=TERMINAL_FIELDS:
        raise AssertionError('Unexpected terminal_bounds.csv columns')
    runtime={'wall_time','wall_seconds'}
    evaluations={'nfev','nfev_refined','sum_nfev'}
    norm_defects={'max_norm_error','norm_defect_refined','max_abs_norm_error','max_norm_defect'}
    parameters={'r','R','Lambda','N','N_gauss','N_kappa','mu','Tmax','N_T','h','L_T','delta_kappa','a','b','K','cutoff_Lambda','Lambda_interval_lo','Lambda_interval_hi','rtol','atol','lambda_','number_of_modes'}
    tight={'F_order0','F_order2','F_refined','interval_integral','leading_tail','tail_leading'}
    checked=0;max_abs=0.
    for row_index,(x,y) in enumerate(zip(a,b),1):
        if terminal:
            terminal_rows.append(compare_terminal_diagnostics(x,y,row_index))
        for key in ah:
            try:u,v=float(x[key]),float(y[key])
            except ValueError:
                if x[key]!=y[key]:raise AssertionError(f'Text differs: {actual.name}:{row_index}:{key}')
                continue
            if not math.isfinite(u) or not math.isfinite(v):
                raise AssertionError(f'Nonfinite result: {actual.name}:{row_index}:{key}')
            if key in runtime:
                if u<0:raise AssertionError(f'Negative elapsed time: {actual.name}')
                continue
            if key in evaluations:
                if u<=0 or not u.is_integer():raise AssertionError(f'Invalid evaluation count: {actual.name}')
                continue
            if key in norm_defects:
                if abs(u)>5e-9:raise AssertionError(f'Unitary norm defect: {actual.name}:{row_index}:{key}: {u}')
                continue
            atol=5e-14 if key in tight else TABLE_ATOL
            ok=((terminal and key in TERMINAL_DERIVED) or
                ((u==v) if key in parameters else math.isclose(u,v,rel_tol=TABLE_RTOL,abs_tol=atol)))
            if not ok:
                raise AssertionError(f'Numerical mismatch {actual.name}:{row_index}:{key}: {u} != {v}')
            checked+=1;max_abs=max(max_abs,abs(u-v))
    result={'rows':len(a),'compared_numeric_fields':checked,'rtol':TABLE_RTOL,'atol':TABLE_ATOL,
            'small_F_and_tail_atol':5e-14,'max_norm_defect':5e-9,
            'performance_fields_excluded':sorted((runtime|evaluations)&set(ah)),
            'norm_defect_fields_bounded':sorted(norm_defects&set(ah)),
            'max_abs_difference':max_abs,'passed':True}
    if terminal:
        result['terminal_bound_diagnostics']={
            'Rs':TERMINAL_RS,'operand_ulps_per_run':TERMINAL_OPERAND_ULPS,
            'ratio_ulps_per_run':TERMINAL_RATIO_ULPS,'bound_relative_slack':TERMINAL_BOUND_SLACK,
            'rows':terminal_rows,
        }
    return result

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
        if not all(math.isfinite(float(r[key])) for key in ('F','norm_error','nfev')) or r['nfev']<=0:
            raise AssertionError(f'Nonfinite or empty transition computation {z}: {r}')
        if abs(r['norm_error'])>5e-9:raise AssertionError(f'Unitary norm defect {z}: {r}')
    if not vals['4.4']['F']<vals['4.6']['F']<vals['5.0']['F']:
        raise AssertionError(f'Nonmonotonicity failed: {vals}')
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=('quick','full'),default='quick')
    ap.add_argument('--report',help='Optional explicit path for JSON report outside the repository')
    a=ap.parse_args()
    report=Path(a.report).expanduser().resolve() if a.report else None
    if report and report.is_relative_to(PAPER):ap.error('--report must be outside the Paper III directory')
    logs=report.parent/(report.stem+'_logs') if report else None
    if logs:logs.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    summary=dict(mode=a.mode,python=sys.version.split()[0],base='Paper III v1.2.0 candidate',checks={},steps=[],passed=False,comments='Scripts are run in a disposable copy, not against the released files.')
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
    env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONUNBUFFERED='1')
    try:
        with tempfile.TemporaryDirectory(prefix='fermions_paperIII_') as t:
            p=Path(t)/'paper_III'
            shutil.copytree(PAPER,p,ignore=shutil.ignore_patterns('__pycache__'))
            summary['checks']=check_frozen_baselines(p)
            if a.mode=='full':
                # Each compared output must be produced during this run. Leaving
                # copied baselines in place could make an omitted write pass.
                for name in GENERATED_TABLES:(p/name).unlink()
            todo=FULL if a.mode=='full' else QUICK
            for i,(label,name) in enumerate(todo,1):
                st=time.perf_counter()
                print(f'[{i}/{len(todo)}] {label}',flush=True)
                path=p/name
                proc=subprocess.run([sys.executable,str(path)],cwd=path.parent,env=env,capture_output=True,text=True)
                elapsed=time.perf_counter()-st
                details={'name':label,'file':name,'returncode':proc.returncode,'elapsed_s':round(elapsed,3),'stdout_tail':proc.stdout[-1800:],'stderr_tail':proc.stderr[-800:]}
                summary['steps'].append(details)
                if logs:
                    (logs/(label+'.out')).write_text(proc.stdout)
                    (logs/(label+'.err')).write_text(proc.stderr)
                if proc.returncode:
                    print(proc.stdout[-3000:],proc.stderr[-3000:],flush=True)
                    raise RuntimeError(f'CHECK FAILED {label}: {proc.returncode}')
                print(f'   PASS ({elapsed:.1f}s)',flush=True)
            if a.mode=='full' and logs:
                generated=logs/'generated_tables'
                summary['generated_tables_dir']=str(generated)
                for name in GENERATED_TABLES:
                    if (p/name).is_file():
                        target=generated/name
                        target.parent.mkdir(parents=True,exist_ok=True)
                        shutil.copyfile(p/name,target)
            st=time.perf_counter()
            summary['fresh_transition_window']=quick_f_check(p)
            print('   Fresh transition window: PASS',round(time.perf_counter()-st,2),'s',flush=True)
            # Recheck newly generated outputs, not only their archived predecessors.
            summary['generated_checks']=check_frozen_baselines(p)
            if a.mode=='full':
                summary['generated_tables']={}
                errors={}
                for name in GENERATED_TABLES:
                    try:
                        summary['generated_tables'][name]=compare_generated_table(p/name,PAPER/name)
                    except (AssertionError,OSError,ValueError,csv.Error) as exc:
                        errors[name]=str(exc)
                if errors:
                    summary['generated_table_errors']=errors
                    raise AssertionError('\n'.join(errors.values()))
            summary['passed']=True
    except Exception as exc:
        summary['error']=f'{type(exc).__name__}: {exc}'
        print(summary['error'],file=sys.stderr,flush=True)
    finally:
        summary['elapsed_s']=round(time.perf_counter()-started,3)
        if report:
            report.parent.mkdir(parents=True,exist_ok=True)
            report.write_text(json.dumps(summary,indent=2)+'\n')
            print('Report:',report,flush=True)
    state='PASS' if summary['passed'] else 'FAIL'
    print(f'Paper III {a.mode.upper()} checks: {state} ({summary["elapsed_s"]:.1f} seconds)',flush=True)
    return 0 if summary['passed'] else 1

if __name__=='__main__':raise SystemExit(main())

#!/usr/bin/env python3
"""Validate the combined candidate in disposable copies. Never stages or publishes.

The original v1.1.0 scientific files are checked by SHA-256 before any execution.
Full mode executes the author's complete independent scalar audit (SKIP_A unset),
recovered AoP refinements, historical fermionic diagnostics, and the 14-stage
fermionic runner. Quick mode is explicitly a subset, not full reproduction.
"""
from __future__ import annotations
from pathlib import Path
import argparse,csv,difflib,hashlib,importlib,json,math,os,re,shutil,subprocess,sys,tempfile,time
from decimal import Decimal,InvalidOperation
ROOT=Path(__file__).resolve().parents[1]

# Exact original package inputs for the two audits eligible for resumption.
# A legacy report did not store these hashes; pinning them here prevents reuse
# after a source/reference edit, even if somebody regenerated SHA256SUMS.
RESUME_AUDIT_HASHES={
    'reproducibility/paper_II/independent_audit/audit_sym.py':'07a34e4d61376cd6839bf74111273bde36085d31c84471b675f6e0b83e68a386',
    'reproducibility/paper_II/independent_audit/audit_num.py':'fb176842ba78c96873661e8e1b0cc8278b8a3e2a42033999790fe4a28eaa9f7d',
    'reproducibility/paper_II/independent_audit/audit_sym_output.txt':'64c8e8532c7696ce8a0f05d6c6998293b4bd6a24add06164c6b3914015b4b733',
    'reproducibility/paper_II/independent_audit/audit_num_output.txt':'ec498473885ad643cf7413bdde1918f1c4b08631c4c9e6fdce582499b5441a6c',
}
LEGACY_STDOUT_ERROR='AssertionError: Independent audit stdout differs beyond elapsed seconds; inspect the saved log'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_manifest(root, manifest, expected):
    """Require complete, duplicate-free coverage, not merely matching listed hashes."""
    seen=set()
    for line in manifest.read_text().splitlines():
        match=re.fullmatch(r'([0-9a-fA-F]{64})  (.+)',line)
        if not match:raise ValueError('Malformed checksum record in '+str(manifest))
        digest,name=match.groups();relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe checksum path: '+name)
        path=(manifest.parent/relative).resolve()
        if not path.is_relative_to(root.resolve()):raise ValueError('Unsafe checksum path: '+name)
        canonical=path.relative_to(root.resolve()).as_posix()
        if canonical in seen:raise ValueError('Duplicate checksum path: '+canonical)
        seen.add(canonical)
        if not path.is_file() or sha(path)!=digest.lower():raise RuntimeError('Checksum mismatch: '+canonical)
    if seen!=set(expected):
        raise RuntimeError(f'Manifest inventory differs: {manifest}: missing={sorted(set(expected)-seen)}, unexpected={sorted(seen-set(expected))}')
    return len(seen)

def validate_integrity(root):
    root=root.resolve()
    paths=json.loads((root/'tools/v120_paths.json').read_text())['paths']
    if not paths or len(paths)!=len(set(paths)):raise ValueError('Empty or duplicate candidate inventory')
    for name in paths:
        relative=Path(name);path=(root/relative).resolve()
        if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(root):
            raise ValueError('Unsafe inventory path: '+name)
        if relative.as_posix()!=name or name.startswith('./'):raise ValueError('Noncanonical inventory path: '+name)
        if not path.is_file():raise FileNotFoundError(name)
    lock=json.loads((root/'tools/v110_lock.json').read_text())
    if not lock['immutable']:raise ValueError('Empty immutable v1.1.0 lock')
    checked=0
    for name,digest in lock['immutable'].items():
        if name not in paths:raise ValueError('Immutable file missing from inventory: '+name)
        if not (root/name).is_file() or sha(root/name)!=digest:
            raise RuntimeError('Immutable v1.1.0 file differs: '+name)
        checked+=1
    counts={'root':check_manifest(root,root/'SHA256SUMS',[n for n in paths if n!='SHA256SUMS'])}
    for folder in ('reproducibility/paper_II','reproducibility/paper_III'):
        manifest=folder+'/SHA256SUMS'
        counts[folder]=check_manifest(root,root/manifest,[n for n in paths if n.startswith(folder+'/') and n!=manifest])
    return {'immutable_v110_files':checked,'root_manifest_files':counts['root'],'manifests':counts,'passed':True}

def normal_runtime(text):return re.sub(r';\s*\d+(?:\.\d+)?s(?=\n|$)','; <elapsed>',text)

def norm_residual_fields(text):
    """Extract only printed B1/B3 norm residuals, preserving every other byte."""
    patterns={
        'B1':(re.compile(r'^  kappa=\s*\d+: '),re.compile(r'(\|alpha\|\^2-\|beta\|\^2-1 = )([^;\s]+)')),
        'B3':(re.compile(r'^  l=\d+ K=[^:\s]+: '),re.compile(r'(norm err )([^;\s]+)')),
    }
    section=None;fields=[];lines=[]
    for line_number,line in enumerate(text.splitlines(keepends=True),1):
        heading=re.match(r'^\(([A-Za-z]\d*)\)',line)
        if heading:section=heading.group(1)
        if section in patterns and patterns[section][0].match(line):
            def replace(match):
                token=match.group(2)
                try:value=Decimal(token)
                except InvalidOperation:raise AssertionError(f'Invalid norm residual at line {line_number}: {token}')
                if not value.is_finite():raise AssertionError(f'Nonfinite norm residual at line {line_number}: {token}')
                fields.append({'section':section,'row':line.split(':',1)[0].strip(),
                               'line':line_number,'printed':token,'value':value})
                return match.group(1)+f'<norm_residual_{len(fields)}>'
            line=patterns[section][1].sub(replace,line)
        lines.append(line)
    return ''.join(lines),fields

def compare_stdout(got,expected,*,numeric_norm_residuals=False,diff_path=None):
    """Exact stdout except elapsed seconds and opt-in B1/B3 norm diagnostics.

    Only the printed norm residual differences have atol=1e-13, rtol=0. The
    printed observables, phases, Wronskians, row IDs and all other text remain
    exact. Decimal arithmetic implements the printed-token boundary exactly.
    """
    a,b=normal_runtime(got),normal_runtime(expected)
    result={'byte_identical':got==expected,'identical_excluding_elapsed_seconds':a==b,
            'numeric_norm_residuals':{'enabled':numeric_norm_residuals,'atol':1e-13 if numeric_norm_residuals else 0.,
                                      'rtol':0.,'fields':[]}}
    try:
        if numeric_norm_residuals:
            masked_a,actual=norm_residual_fields(a);masked_b,reference=norm_residual_fields(b)
            if masked_a!=masked_b or len(actual)!=len(reference):
                raise AssertionError('Text or scientific values differ outside the B1/B3 norm-residual fields')
            for x,y in zip(actual,reference):
                delta=x['value']-y['value'];absolute=abs(delta)
                record={'section':y['section'],'row':y['row'],'line':y['line'],
                        'reference_printed':y['printed'],'actual_printed':x['printed'],
                        'reference':float(y['value']),'actual':float(x['value']),
                        'delta':float(delta),'abs_delta':float(absolute)}
                result['numeric_norm_residuals']['fields'].append(record)
                if absolute>Decimal('1e-13'):
                    raise AssertionError(f'Norm-residual difference exceeds 1e-13 at {y["section"]} {y["row"]}: reference={y["printed"]}, actual={x["printed"]}, delta={delta}')
            result['numeric_norm_residuals']['max_abs_delta']=max((x['abs_delta'] for x in result['numeric_norm_residuals']['fields']),default=0.)
        elif a!=b:
            raise AssertionError('Text differs beyond elapsed seconds')
    except AssertionError as exc:
        location=''
        if diff_path is not None:
            diff_path=Path(diff_path)
            diff_path.parent.mkdir(parents=True,exist_ok=True)
            header=(f'Reference length: {len(b)}; actual length: {len(a)}\n'
                    f'Reference final newline: {b.endswith(chr(10))}; actual final newline: {a.endswith(chr(10))}\n'
                    f'Numeric norm policy enabled: {numeric_norm_residuals}; atol=1e-13, rtol=0\n'
                    f'Reason: {exc}\n\n')
            diff=''.join(difflib.unified_diff(b.splitlines(keepends=True),a.splitlines(keepends=True),fromfile='reference (elapsed normalized)',tofile='actual (elapsed normalized)'))
            diff_path.write_text(header+diff)
            location=f'; saved diff: {diff_path}'
        raise AssertionError(f'Independent audit stdout comparison failed: {exc}{location}') from exc
    result['passed']=True
    return result

def validate_resume_destinations(source,dest,logs):
    """A resumed validation must preserve its original report and logs."""
    source=source.expanduser().resolve();dest=dest.resolve();logs=logs.resolve()
    original_logs=(source.parent/(source.stem+'_logs')).resolve()
    if dest==source or dest.is_relative_to(original_logs):
        raise ValueError('Resume requires a new report outside the original report/logs')
    if logs==original_logs or logs.is_relative_to(original_logs) or original_logs.is_relative_to(logs):
        raise ValueError('Resume requires a separate new log directory')
    if dest.exists() or logs.exists():
        raise ValueError('Resume destination report/logs already exist; choose a fresh report path')

def load_resume(source,*,mode,paper,environment,root=ROOT):
    """Validate the narrowly allowed legacy failure and load only its two audits."""
    source=source.expanduser().resolve();root=root.resolve()
    old=json.loads(source.read_text());old_logs=source.parent/(source.stem+'_logs')
    if mode!='full' or paper not in ('II','all') or old.get('mode')!=mode or old.get('paper')!=paper:
        raise ValueError('Resume requires identical full mode and Paper II/all selection')
    if old.get('passed') is not False or old.get('error')!=LEGACY_STDOUT_ERROR:
        raise ValueError('Resume is allowed only for the original exact-stdout comparison failure')
    if old.get('integrity',{}).get('passed') is not True:
        raise ValueError('Original report does not record successful input integrity')
    if old.get('environment')!=environment:
        raise ValueError('Resume environment differs from the original report')
    names=['paperII_independent_sym','paperII_independent_num']
    steps=old.get('steps',[])
    if [s.get('name') for s in steps]!=names or any(type(s.get('returncode')) is not int or s['returncode']!=0 for s in steps):
        raise ValueError('Original report must contain exactly the two successful independent audit executions')
    comparisons=old.get('comparisons',{})
    if set(comparisons)!={'audit_sym'} or comparisons['audit_sym'].get('identical_excluding_elapsed_seconds') is not True:
        raise ValueError('Original report does not record the expected symbolic comparison')
    for relative,digest in RESUME_AUDIT_HASHES.items():
        path=root/relative
        if not path.is_file() or sha(path)!=digest:
            raise ValueError('Resume source/reference differs from the original package: '+relative)
    cache={};records={}
    for name,short in zip(names,('sym','num')):
        out=old_logs/(name+'.out');err=old_logs/(name+'.err')
        if not out.is_file() or not err.is_file():
            raise FileNotFoundError('Missing original stdout/stderr for '+name)
        relative=f'reproducibility/paper_II/independent_audit/audit_{short}.py'
        ref=f'reproducibility/paper_II/independent_audit/audit_{short}_output.txt'
        proof={'reused_from':str(source),'original_stdout':str(out),'original_stderr':str(err),
               'stdout_sha256':sha(out),'stderr_sha256':sha(err),
               'source':relative,'source_sha256':sha(root/relative),
               'reference':ref,'reference_sha256':sha(root/ref)}
        cache[name]={'stdout':out.read_text(),'stderr':err.read_text(),'provenance':proof}
        records[name]=proof
    meta={'reused_from':str(source),'original_report_sha256':sha(source),
          'allowed_cached_steps':names,'artifacts':records,
          'policy':'Original full-mode failure and environment verified; original package source/reference hashes pinned; saved stdout is compared again. No other step is reused.'}
    return cache,meta

JOST_IR_NORM_ERROR_LIMIT = Decimal('1e-6')


def csv_compare(actual, expected, *, jost_ir=False):
    """Compare AoP tables with a bounded norm defect in the Jost IR diagnostic.

    The Jost coefficients grow as O(K**(-1/2)), and their squared magnitudes
    cancel to a unit Wronskian. Bound its absolute defect by one part per million
    and compare the scaled amplitude with the standard CSV tolerances.
    """
    if jost_ir and (actual.name != 'Jost_IR_diagnostic.csv' or expected.name != actual.name):
        raise AssertionError('Jost IR norm policy requires two Jost_IR_diagnostic.csv files')
    with actual.open(newline='') as f:
        reader=csv.DictReader(f);a=list(reader);aheader=reader.fieldnames
    with expected.open(newline='') as f:
        reader=csv.DictReader(f);b=list(reader);bheader=reader.fieldnames
    if not aheader or aheader!=bheader or len(aheader)!=len(set(aheader)):
        raise AssertionError(f'Columns differ or are empty/duplicated: {actual.name}')
    if jost_ir and aheader != ['K', 'sqrtK_alpha_plus_beta', 'norm_error']:
        raise AssertionError(f'Unexpected Jost IR columns: {aheader}')
    if not a or not b:raise AssertionError(f'Empty results: {actual.name}')
    if len(a)!=len(b):raise AssertionError(f'Row count differs: {actual.name}')
    max_abs=0.;n=0
    norm_residuals = []
    for x,y in zip(a,b):
        if None in x or None in y or any(v is None for v in [*x.values(),*y.values()]):
            raise AssertionError(f'Malformed CSV row: {actual.name}')
        if x.keys()!=y.keys():raise AssertionError(f'Columns differ: {actual.name}')
        for k in x:
            try:u,v=float(x[k]),float(y[k])
            except ValueError:
                if jost_ir:
                    raise AssertionError(f'Invalid numeric result {actual.name}:{k}: {x[k]}, {y[k]}')
                if x[k]!=y[k]:raise AssertionError(f'Text differs: {actual.name}:{k}')
                continue
            if not math.isfinite(u) or not math.isfinite(v):
                raise AssertionError(f'Nonfinite result {actual.name}:{k}: {u}, {v}')
            if jost_ir and k == 'norm_error':
                if any(not 0 <= Decimal(token) <= JOST_IR_NORM_ERROR_LIMIT for token in (x[k], y[k])):
                    raise AssertionError(
                        f'{actual.name}:norm_error at K={x["K"]} outside '
                        f'[0, {JOST_IR_NORM_ERROR_LIMIT}] (actual {u}, reference {v})'
                    )
                norm_residuals.append({'K': float(x['K']), 'actual': u, 'reference': v,
                                       'abs_delta': abs(u-v)})
                continue
            if jost_ir and k == 'K':
                if Decimal(x[k]) != Decimal(y[k]):
                    raise AssertionError(f'Frequency mismatch {actual.name}:K: {x[k]} != {y[k]}')
            elif not math.isclose(u,v,rel_tol=5e-9,abs_tol=5e-11):
                raise AssertionError(f'Numerical mismatch {actual.name}:{k}: {u} != {v}')
            max_abs=max(max_abs,abs(u-v));n+=1
    result = {'rows':len(a),'numeric_fields':n,'max_abs_difference':max_abs,'rtol':5e-9,'atol':5e-11}
    if jost_ir:
        result['jost_ir_norm_error'] = {
            'absolute_bound': float(JOST_IR_NORM_ERROR_LIMIT),
            'max_actual': max(row['actual'] for row in norm_residuals),
            'max_reference': max(row['reference'] for row in norm_residuals),
            'rows': norm_residuals,
        }
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=['quick','full'],default='quick')
    ap.add_argument('--paper',choices=['II','III','all'],default='all')
    ap.add_argument('--report',type=Path,required=True,help='Report outside the repository; logs are saved alongside it')
    ap.add_argument('--resume-from',type=Path,help='Original failed full report; reuse only its two verified independent-audit stdout logs')
    args=ap.parse_args();dest=args.report.expanduser().resolve()
    if dest.is_relative_to(ROOT.resolve()):ap.error('--report must be outside the repository')
    logs=dest.parent/(dest.stem+'_logs')
    if args.resume_from:
        try:validate_resume_destinations(args.resume_from,dest,logs)
        except ValueError as exc:ap.error(str(exc))
    dest.parent.mkdir(parents=True,exist_ok=True);logs.mkdir(exist_ok=True)
    report={'mode':args.mode,'paper':args.paper,'published':False,'steps':[],'comparisons':{},'passed':False}
    started=time.monotonic()
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONUNBUFFERED='1');env.pop('SKIP_A',None)
    # A caller's PYTHONOPTIMIZE must never disable scientific assertions.
    env.pop('PYTHONOPTIMIZE',None)
    resume_cache={}
    def run(name,path,extra=()):
        if name in resume_cache:
            item=resume_cache.pop(name)
            if sha(path)!=item['provenance']['source_sha256']:
                raise RuntimeError('Cached audit source does not match the verified copy: '+name)
            print('REUSE',name,'from',item['provenance']['original_stdout'],flush=True)
            (logs/(name+'.out')).write_text(item['stdout']);(logs/(name+'.err')).write_text(item['stderr'])
            report['steps'].append({'name':name,'returncode':0,'elapsed_s':0.,'execution':'reused_saved_output',**item['provenance']})
            print('DONE',name,'(code=0; reused output; comparison pending)',flush=True)
            return item['stdout']
        print('RUN',name,flush=True);st=time.monotonic()
        proc=subprocess.run([sys.executable,str(path),*map(str,extra)],cwd=path.parent,env=env,capture_output=True,text=True)
        (logs/(name+'.out')).write_text(proc.stdout);(logs/(name+'.err')).write_text(proc.stderr)
        report['steps'].append({'name':name,'returncode':proc.returncode,'elapsed_s':time.monotonic()-st,
                                'execution':'fresh','source_sha256':sha(path),
                                'stdout_sha256':sha(logs/(name+'.out')),'stderr_sha256':sha(logs/(name+'.err'))})
        if proc.returncode:raise RuntimeError(f'{name} failed; see {logs}')
        print('DONE',name,'(code=0; validation pending)',flush=True);return proc.stdout
    try:
        report['integrity']=validate_integrity(ROOT)
        report['environment']={'python':sys.version,**{n:importlib.import_module(n).__version__ for n in ['numpy','scipy','pandas','sympy','mpmath']}}
        if args.resume_from:
            resume_cache,report['resume']=load_resume(args.resume_from,mode=args.mode,paper=args.paper,environment=report['environment'],root=ROOT)
        with tempfile.TemporaryDirectory(prefix='schwarzschild_v120_') as td:
            temp=Path(td)
            if args.paper in ('II','all'):
                p=temp/'paper_II';shutil.copytree(ROOT/'reproducibility/paper_II',p,ignore=shutil.ignore_patterns('__pycache__','results'))
                a=p/'independent_audit';s=p/'aop_supplement';(s/'results').mkdir(exist_ok=True)
                got=run('paperII_independent_sym',a/'audit_sym.py')
                report['comparisons']['audit_sym']=compare_stdout(got,(a/'audit_sym_output.txt').read_text(),diff_path=logs/'paperII_independent_sym.diff.txt')
                print('PASS paperII_independent_sym (stdout compared)',flush=True)
                if args.mode=='full':
                    got=run('paperII_independent_num',a/'audit_num.py')
                    report['comparisons']['audit_num']=compare_stdout(got,(a/'audit_num_output.txt').read_text(),numeric_norm_residuals=True,diff_path=logs/'paperII_independent_num.diff.txt')
                    print('PASS paperII_independent_num (stdout and B1/B3 norm residuals compared)',flush=True)
                    run('paperII_AoP_revision',s/'check_AoP_revision.py')
                    run('paperII_publication_tables',s/'export_publication_tables.py')
                    aop_snapshot = logs / 'paperII_aop_results'
                    shutil.copytree(s/'results', aop_snapshot, dirs_exist_ok=True)
                    report['aop_results_dir'] = str(aop_snapshot)
                    comparison_errors = {}
                    for name in ['swave_refinement.csv','normalized_packet.csv','packet_refinement.csv','Jost_IR_diagnostic.csv','squeezed_pair.csv','table1_swave_AoP.csv','table2_AoP.csv']:
                        try:
                            report['comparisons'][name] = csv_compare(
                                s/'results'/name, s/'reference_results'/name,
                                jost_ir=(name == 'Jost_IR_diagnostic.csv'),
                            )
                        except AssertionError as exc:
                            comparison_errors[name] = str(exc)
                    if comparison_errors:
                        report['aop_comparison_errors'] = comparison_errors
                        raise AssertionError(
                            '\n'.join(comparison_errors.values())
                            + f'\nFresh AoP tables saved in {aop_snapshot}'
                        )
                    aj=json.loads((s/'results/audit_summary.json').read_text())
                    if aj.get('count')!=29 or not aj.get('all_pass'):raise AssertionError('AoP assertions failed')
                    report['aop_29_checks']=aj
                    print('PASS paperII_AoP_revision and publication tables (29 assertions and seven CSV comparisons)',flush=True)
                run('paperII_scalar_Ward',s/'check_scalar_Ward.py')
                wj=json.loads((s/'results/scalar_Ward_jets.json').read_text())
                if not wj.get('all_pass') or len(wj['cases'])!=20:raise AssertionError('Ward checks failed')
                report['scalar_Ward_20']=wj
                print('PASS paperII_scalar_Ward (20 cases)',flush=True)
            if args.paper in ('III','all'):
                p=temp/'paper_III';shutil.copytree(ROOT/'reproducibility/paper_III',p,ignore=shutil.ignore_patterns('__pycache__'))
                internal=logs/'paperIII_runner.json'
                run('paperIII_'+args.mode,p/'run_checks.py',['--mode',args.mode,'--report',internal])
                j=json.loads(internal.read_text());report['paperIII']=j
                expected_steps=14 if args.mode=='full' else 5
                if j.get('passed') is not True or j.get('mode')!=args.mode or len(j.get('steps',[]))!=expected_steps:
                    raise AssertionError('Paper III runner failed or did not execute every requested step')
                print('PASS paperIII_'+args.mode+' (runner report checked)',flush=True)
                if args.mode=='full':
                    l=p/'legacy_phase2'
                    for name in ['dirac_phase0.py','dirac_phase2_numerics.py','dirac_phase2_packet.py','dirac_stress_identity_check.py']:
                        run('legacy_'+name[:-3],l/name)
                    # These archived scripts print diagnostics rather than asserting all values.
                    report['legacy_note']='Four scripts executed; historical auxiliary sign/finite-mode scope documented in their README. Original provided outputs are retained.'
        if resume_cache:raise RuntimeError('A permitted cached audit was not consumed')
        report['integrity_after']=validate_integrity(ROOT)
        report['passed']=True
    except Exception as exc:
        report['error']=f'{type(exc).__name__}: {exc}'
        print(report['error'],file=sys.stderr)
    finally:
        report['elapsed_s']=time.monotonic()-started
        dest.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        print('Report:',dest,flush=True)
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())

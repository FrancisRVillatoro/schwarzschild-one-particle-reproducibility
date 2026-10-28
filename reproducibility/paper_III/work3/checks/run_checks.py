#!/usr/bin/env python3
"""Run the Work3 checks with the active Python interpreter; stop on any failure."""
import json,os,platform,subprocess,sys
from pathlib import Path
import numpy,scipy,sympy
root=Path(__file__).resolve().parent
env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1';env['OPENBLAS_NUM_THREADS']='1'
names=['check_work3_symbolics.py','check_invariant_harmonics.py','check_continuum_stress.py','check_continuum_refinements.py']
for i,name in enumerate(names,1):
 print(f'[{i}/{len(names)}] {name}',flush=True)
 completed=subprocess.run([sys.executable,str(root/name)],cwd=root,env=env,check=True,text=True,capture_output=True)
 (root/(Path(name).stem+'_output.txt')).write_text(completed.stdout)
 print(completed.stdout,flush=True)
versions={'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__,'tests':names,'all_process_return_codes':0}
(root/'environment_and_run.json').write_text(json.dumps(versions,indent=2)+'\n')
print('All Work3 checks completed successfully.',flush=True)

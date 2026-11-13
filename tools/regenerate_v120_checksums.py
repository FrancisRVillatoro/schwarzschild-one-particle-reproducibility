#!/usr/bin/env python3
"""Regenerate checksums for the explicit candidate inventory, never a recursive .venv scan."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    lock=json.loads((ROOT/'tools/v110_lock.json').read_text())
    for n,h in lock['immutable'].items():
        if digest(ROOT/n)!=h:raise RuntimeError('Immutable v1.1.0 file changed: '+n)
    paths=json.loads((ROOT/'tools/v120_paths.json').read_text())['paths']
    for name in paths:
        p=(ROOT/name).resolve()
        if not p.is_relative_to(ROOT.resolve()):raise ValueError('Unsafe inventory path')
        if name.endswith('/SHA256SUMS') or name=='SHA256SUMS':continue
        if not p.is_file():raise FileNotFoundError(name)
    for folder in ['reproducibility/paper_II','reproducibility/paper_III']:
        records=[n for n in paths if n.startswith(folder+'/') and n!=folder+'/SHA256SUMS']
        text=''.join(f'{digest(ROOT/n)}  {n[len(folder)+1:]}\n' for n in sorted(records))
        (ROOT/folder/'SHA256SUMS').write_text(text)
    text=''.join(f'{digest(ROOT/n)}  ./{n}\n' for n in sorted(paths) if n!='SHA256SUMS')
    (ROOT/'SHA256SUMS').write_text(text)
    print('Regenerated root, Paper-II and Paper-III manifests; original scientific files unchanged.')
if __name__=='__main__':main()

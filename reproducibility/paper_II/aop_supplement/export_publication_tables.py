#!/usr/bin/env python3
"""Export the documented refined rows; never alter v1.1.0 production scripts.

Table 1 top uses r_min=1e-14, rtol=2e-13, atol=2e-15.
Table 2 uses 96 Gauss-Legendre nodes and the same fine tolerance.
The lower half of Table 1 remains the fixed-tower run_part1.py output,
independently reproduced by independent_audit/audit_num.py.
"""
from pathlib import Path
import argparse,csv,json,hashlib

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def write(path, rows):
    if not rows:raise ValueError(f'No selected rows for {path.name}')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,default=Path(__file__).resolve().parent/'results');a=ap.parse_args();d=a.results
    sw=[r for r in read(d/'swave_refinement.csv') if float(r['rmin'])==1e-14 and float(r['rtol'])==2e-13]
    if len(sw)!=7:raise ValueError('Expected seven refined s-wave entries')
    write(d/'table1_swave_AoP.csv',[dict(K=int(float(r['K'])),K_beta2=r['K_beta2'],printed_6dp=f"{float(r['K_beta2']):.6f}",rmin=r['rmin'],rtol=r['rtol'],atol=2e-15) for r in sw])
    pk=[r for r in read(d/'normalized_packet.csv') if int(r['N_omega'])==96 and float(r['rtol'])==2e-13]
    out=[]
    for radius in (1e-2,1e-4,1e-6,1e-8):
        row={'r':radius}
        for loc in ('Y00','pole','equator'):
            for energy in (0.,.5,2.):
                matches=[r for r in pk if float(r['r'])==radius and r['location']==loc and float(r['E'])==energy]
                if len(matches)!=1:raise ValueError((radius,loc,energy,len(matches)))
                r=matches[0];key=f'{loc}_E{energy:g}'
                val=float(r['r3rho'] if loc=='equator' else r['ratio'])
                row[key]=val
                row[key+'_printed']=f'{val:.3e}' if loc=='equator' else f'{val:.6f}'
        out.append(row)
    write(d/'table2_AoP.csv',out)
    data={'Table_1_top':{'source':'swave_refinement.csv','rmin':1e-14,'rtol':2e-13,'atol':2e-15},'Table_1_bottom':{'source':'../run_part1.py','independent_check':'../independent_audit/audit_num.py','unchanged_production':True},'Table_2':{'source':'normalized_packet.csv','N_omega':96,'rtol':2e-13,'atol':2e-15,'Z':float(pk[0]['Z']),'scope':'Normalized horizon-flux diagnostic, not exterior Unruh scattering'},'source_sha256':{n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in ['swave_refinement.csv','normalized_packet.csv','packet_refinement.csv']}}
    (d/'publication_table_sources.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: exported Table 1 top and Table 2; Z=',data['Table_2']['Z'])
if __name__=='__main__':main()
